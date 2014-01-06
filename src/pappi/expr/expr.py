from .. import sql
from ..table_manager import TableManager
from ..utils import progressbar
import csv

# at least 90 percent of genes have to be covered
# TODO:
# (i should check what the optimal thresholds are)
PARAM_TISSUE_MIN_GENE_COVERAGE = 90


class GeneExpression(TableManager):

    file_field_seperator = '\t'
    file_has_header = True
    file_quoting = csv.QUOTE_NONE

    """
    An interface/super-class to all Gene/Protein-Expression data sets.
    """
    def __init__(self, filename, sql_connection, dataset_name):
        """
        Constructs a new GeneExpression class. Should not be used
        directly but only from overwriting subclasses.
        """
        self.filename = filename
        self.sql_conn = sql_connection
        self.name = dataset_name

    def init_data(self):
        # TODO re-evaluate which order is best ... !?
        self.import_raw_file()
        self.linearize_table()
        # filter before ID mapping, because if two (or more) IDs are mapped to
        # the same ID and then combined, only `reliable`/`filtered` rows
        # should be used
        self.filter()
        self.id_mapping()
        self.normalize_table()
        # step to group and aggregate duplicate rows
        self.rm_duplicates()
        self.classify()

        # summarize the genes/proteins with stats of how much
        # they are expressed (expressed_count / total_count)
        self.expr_counts()

    def import_raw_file(self):
        """
        Imports the expression file in raw format into the SQL database
        If the specific datasets need some kind of raw pre-processing,
        then this function MUST be overwritten.
        """
        table_name = self.next_tmp_table("raw")
        sql.import_csv(self.filename, table_name, self.file_field_seperator,
                       self.file_has_header, csv_quoting=self.file_quoting,
                       sql_conn=self.sql_conn)

    def linearize_table(self):
        """
        In case the table structure is [Gene] x [Tissue/Cell-Type], the table
        needs to be linearized to a column structure
        [Gene],[Tissue/Cell-Type],[Expr]
        """
        pass

    def filter(self):
        """
        Filters the expression data for genes with certain
        reliablilities or confidences.
        """
        pass

    def id_mapping(self):
        """
        Maps the gene/protein identifier to a unifying identifier.
        """
        pass

    def normalize_table(self):
        """
        Normalizes the table structure to the format
        [Gene],[Type] (,[Specific Type],...), [ExpressionValue]
        """
        pass

    def rm_duplicates(self):
        """
        Removes duplicate rows (caused by ID mapping) and aggregates their
        expression values via MAX(). The maximum expression value is choosen,
        because the expression level of the Gene (from all it's mapped rows) is
        dominated by this maximum.
        """
        # NOTE: this implementation works only for a numerical
        # `ExpressionValue`, in case this column is something other than
        # numeric (i.e. a string), this function has to be overwritten by the
        # subclasses
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table("no_dupl")
        # TODO: count how many are removed, and log somewhere
        sqlquery = ('SELECT Gene, Type, '
                    'MAX(ExpressionValue) AS ExpressionValue '
                    'FROM ' + src_table + ' '
                    'GROUP BY Gene, Type')
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def classify(self):
        """
        Thresholds or classifies genes into either expressed ( = 1)
        or non-expressed (= 0).
        """
        if not hasattr(self, 'classify_cond'):
            raise NotImplementedError("The classify method has to be "
                                      "implemented by all subclasses, or "
                                      "the self.classify_cond attribute set")
        # get the src and dest table, assuming they are in normalized format,
        # create the expression-classified version of the table
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table("")
        sqlquery = ('SELECT Gene, Type, '
                    ' CASE WHEN ExpressionValue ' + self.classify_cond + ' '
                    ' THEN 1 ELSE 0 END AS Expressed '
                    'FROM ' + src_table)
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def expr_counts(self):
        """
        Creates a table with the total counts and expressed counts of each
        gene/protein in the expression data set. This is a precursor for
        tissue-specific v.s. housekeeping classifications.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table("expr_counts")
        sqlquery = ('SELECT Gene, '
                    ' COUNT(Type) AS TotalCount, '
                    ' SUM(Expressed) AS ExpressedCount '
                    'FROM ' + src_table)
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def create_ids_table(self):
        """
        Creates a table named <name>_ids that holds all the distinct IDs used
        in the expression data set, where <name> is the name of the expression
        data set.
        """
        if sql.table_exists(self.name + '_ids', self.sql_conn):
            return
        cur = sql_conn.cursor()
        cur.execute('CREATE TABLE IF NOT EXISTS `' + self.name + '_ids` '
                    '(id integer primary key autoincrement, Gene varchar(16))')
        sqlquery = ('SELECT DISTINCT Gene FROM ' + self.name)
        cur.execute('INSERT INTO `' + self.name + '_ids` (Gene) ' + sqlquery)
        cur.close()
        sql_conn.commit()

    def create_tissue_table(self):
        """
        Creates a table <name>_tissues with all unique tissue/cell types in the
        expression data set in sorted order, where <name> is the name of the
        expression data set.
        """
        if sql.table_exists(self.name + '_tissues', self.sql_conn):
            return
        sqlquery = ('SELECT DISTINCT Type, '
                    'COUNT(DISTINCT Gene) AS Gene_Coverage '
                    'FROM ' + self.name + ' GROUP BY Type ORDER BY Type')
        # sqlquery = ('SELECT DISTINCT Type FROM ' + self.name + ' '
        #            'ORDER BY Type')
        sql.new_table_from_query(self.name + '_tissues', sqlquery,
                                 self.sql_conn)

    def create_tissue_coverage_table(self):
        """
        """

        # create tissue table first
        self.create_tissue_table()

        # get SQL cursor
        cur = self.sql_conn.cursor()

        # first create the results table
        cur.execute('DROP TABLE IF EXISTS ' + self.name + '_coverage')
        cur.execute('CREATE TABLE ' + self.name + '_coverage '
                    '(gene_coverage_threshold int, gene_coverage int, '
                    'total_genes int, tissue_coverage int, total_tissues int)')

        # get the tissue table with the gene coverage for each tissue
        cur.execute('SELECT Type, Gene_Coverage FROM ' + self.name +
                    '_tissues ORDER BY Gene_Coverage ASC')
        all_rows = cur.fetchall()
        num_tissues = len(all_rows)
        cur.execute('SELECT COUNT(DISTINCT Gene) FROM ' + self.name)
        num_genes_total = cur.fetchone()[0]
        print("Creating coverage table:")
        progressbar.show_progress(0.0)

        # TODO:
        #       formalize this `greedy` algorithm. There might be a more
        #       optimal combination of tissues and proteins that cover more
        #       in total. I'm guessing this problem is NP hard (not sure
        #       though, i may have to check that out) also the greedy algo
        #       could be extended, so that it adds further tissues even if the
        #       next one doesn't improve the score (skip one), that would
        #       require though a different more dynamic approach/implementation
        #       of the greedy algorithm, and most possibly not via SQL but with
        #       python using 2D matrizzes instead of a ~ 1D table for Genes
        #       CROSS Tissues
        for coverage in sorted(set([int(c) for t, c in all_rows])):
            num_bigger = sum(map(lambda x: x[1] >= coverage, all_rows))
            sqlquery = ('SELECT COUNT(*) FROM ('
                        'SELECT a.Gene, COUNT(*) FROM ' + self.name + ' AS a '
                        'INNER JOIN ' + self.name + '_tissues AS b ON '
                        'a.Type = b.Type WHERE Gene_Coverage >= '
                        + str(coverage) + ' GROUP BY Gene HAVING COUNT(*) >= '
                        + str(num_bigger)
                        + ')')
            cur.execute(sqlquery)
            num_genes_covered = cur.fetchone()[0]

            cur.execute('INSERT INTO ' + self.name + '_coverage '
                        '(gene_coverage_threshold, gene_coverage, '
                        'total_genes, tissue_coverage, '
                        'total_tissues) '
                        ' VALUES (?,?,?,?,?)',
                        [coverage, num_genes_covered, num_genes_total,
                         num_bigger, num_tissues])
            progressbar.show_progress((num_tissues - num_bigger)
                                      * 1.0 / num_tissues)

           # print("with %.0f coverage threshold: %.2f tissue and %.2f genes"
           #       "are left with full coverage" %
           #       (coverage, num_bigger * 100.0 / num_tissues,
           #        num_genes_covered * 100.0 / num_genes_total))
           # print("num bigger: %i " % num_bigger)
           # print("Coverage is: %.2f" % (coverage))
        progressbar.finish_progress()

    def get_optimal_coverage_threshold(self):
        """
        """
        # TODO this is not _optimal_ per se, an optimal solution is probably
        # NP-hard and very related to the Netflix Price Challenge

        # create the coverage table if it does not yet exists
        # (this may take a while)
        if not sql.table_exists(self.name + '_coverage', self.sql_conn):
            self.create_tissue_coverage_table()

        # get the threshold for maximum coverage
        sqlquery = ('SELECT MIN(gene_coverage_threshold) '
                    'FROM ' + self.name + '_coverage '
                    'WHERE gene_coverage*tissue_coverage = '
                    '(SELECT MAX(gene_coverage*tissue_coverage) '
                    'FROM ' + self.name + '_coverage)')
        cur = self.sql_conn.cursor()
        cur.execute(sqlquery)
        threshold = cur.fetchone()[0]
        # return the result
        return threshold

    def create_core_table(self, threshold=None):
        """
        """
        # if threshold is not give, use "optimal" value
        if threshold is None:
            threshold = self.get_optimal_coverage_threshold()

        # join/intersect the main table with list of covered genes and
        # list of covered tissues

        # get number of tissues in core
        cur = self.sql_conn.cursor()
        cur.execute('SELECT COUNT(*) FROM ' + self.name + '_tissues '
                    'WHERE Gene_Coverage >= ' + str(threshold))
        num_tissues = cur.fetchone()[0]

        # join table with tissues, then delete Genes that are not in the core
        sqlquery = ('SELECT Gene, Type, Expressed FROM ' + self.name + ' '
                    'WHERE Type IN (SELECT Type FROM ' + self.name + '_tissues'
                    ' WHERE Gene_Coverage >= ' + str(threshold) + ')')
        sql.new_table_from_query(self.name + '_core', sqlquery, self.sql_conn)
        print("num tissues: " + str(num_tissues))
        cur.execute('DELETE FROM ' + self.name + '_core '
                    'WHERE Gene IN ('
                    '  SELECT Gene FROM ' + self.name + '_core '
                    '  GROUP BY Gene '
                    '  HAVING COUNT(*) < ' + str(num_tissues) + ''
                    ')')

        cur.close()
        self.sql_conn.commit()

    def get_tissue_table(self, min_coverage=None):
        """
        Returns a SQL query which in turn returns all the tissues that have
        data for at least the `min_coverage` percentage of total genes in the
        data set.

        @param min_coverage:    Only those tissues are returned, that cover at
                                least this percentage of Genes in the
                                expression data set. If this is set to `None`,
                                then ALL tissues are returned. (default: None)
        """
        # first create the tissue table if it isn't yet created
        self.create_tissue_table()

        if min_coverage is None:
            # simply return all Tissues
            return self.name + '_tissues'
        else:
            # return only those that have minimal coverage of genes
            sqlquery = ('(SELECT Type FROM ' + self.name + '_tissues '
                        'WHERE Gene_Coverage >= ' + str(min_coverage) + ')')
            return sqlquery

    def create_node_labels(self, only_complete=True, sep='|', null_syb='-'):
        """
        Creates a table of labels for the expression data set. For each gene
        this creates a label of the expression value for each tissue {0,1}
        concatenated with `sep` as separator and `null_syb` as replacement
        of NULL values.

        One such label has the format: 1|1|0|1|1|1|0|1|-|1|-|0|1|1|0|...
        with one 'column' {0,1,-} for each tissue in the same order as in the
        `name`_tissues table.

        @param only_complete:   Export only labels where there is data
                                available for ALL tissue types (i.e. NULL
                                values are not exported).
        @param sep:             The separator inserted between concatenated
                                values. (default: '|')
        @param null_syb:        The character inserted where no expression
                                value is available (i.e. NULL). default: '-'.
        """
        # first of all, create the tissue table (if it doesn't yet exists)
        tissue_tbl = self.get_tissue_table(PARAM_TISSUE_MIN_GENE_COVERAGE)
        # then create the ids table, if it does not yet exist
        self.create_ids_table()

        sqlquery = ('SELECT a.Gene AS Gene, '
                    'group_concat(CASE WHEN b.Expressed IS NULL '
                    'THEN "' + null_syb + '" ELSE b.Expressed END, '
                    '"' + sep + '") AS Label '
                    'FROM '
                    '(SELECT a.Gene, b.Type FROM ' + self.name + '_ids AS a, '
                    ' ' + tissue_tbl + ' AS b) AS a '
                    ' LEFT OUTER JOIN ' + self.name + ' AS b '
                    ' ON a.Gene = b.Gene AND a.Type = b.Type GROUP BY a.Gene')
        sql.new_table_from_query(self.name + '_node_labels', sqlquery,
                                 self.sql_conn)

    def export_node_labels(self, node_ids_tbl, filename,
                           sep='|', null_syb='-'):
        """
        Exports binary expression node labels for the given Node-IDs table
        which must be of the form (id, Gene).
        """
        # create node labels if they don't yet exist
        if not sql.table_exists(self.name + '_node_labels', self.sql_conn):
            self.create_node_labels(sep, null_syb)

        # map the node labels to the node ids of the given table
        cur = self.sql_conn.cursor()
        cur.execute('SELECT b.id, a.Label '
                    'FROM ' + self.name + '_node_labels AS a '
                    'INNER JOIN ' + node_ids_tbl + ' AS b ON a.Gene = b.Gene')

        # save results to file
        with open(filename, "w") as f:
            for row in cur.fetchall():
                f.write("%i\t%s\n" % (row[0], row[1]))
