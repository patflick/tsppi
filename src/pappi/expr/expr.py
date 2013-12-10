from .. import sql
from ..table_manager import TableManager
import csv


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
        self.id_mapping()
        self.filter()
        self.normalize_table()
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

    def normalize_table(self):
        """
        Normalizes the table structure to the format
        [Gene],[Type] (,[Specific Type],...), [ExpressionValue]
        and especially gets rid of the [Gene] x [Type] format.
        """
        pass

    def id_mapping(self):
        """
        Maps the gene/protein identifier to a unifying identifier.
        """
        pass

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
        Creates a table named `expr`_ids that holds all the distinct IDs used
        in the expression data set.
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
        Creates a table `name`_tissues with all unique tissue/cell types in the
        expression data set in sorted order..
        """
        if sql.table_exists(self.name + '_tissues', self.sql_conn):
            return
        sqlquery = ('SELECT DISTINCT Type FROM ' + self.name + ' '
                    'ORDER BY Type')
        sql.new_table_from_query(self.name + '_tissues', sqlquery,
                                 self.sql_conn)

    def create_node_labels(self, sep='|', null_syb='-'):
        """
        Creates a table of labels for the expression data set. For each gene
        this creates a label of the expression value for each tissue {0,1}
        concatenated with `sep` as separator and `null_syb` as replacement
        of NULL values.

        One such label has the format: 1|1|0|1|1|1|0|1|-|1|-|0|1|1|0|...
        with one 'column' {0,1,-} for each tissue in the same order as in the
        `name`_tissues table.

        @param sep:         The separator inserted between concatenated values.
                            default: '|'.
        @param null_syb:    The character inserted where no expression value
                            is available (i.e. NULL). default: '-'.
        """
        # first of all, create the tissue table (if it doesn't yet exists)
        self.create_tissue_table()
        # then create the ids table, if it does not yet exist
        self.create_ids_table()

        sqlquery = ('SELECT a.Gene AS Gene, '
                    'group_concat(CASE WHEN b.Expressed IS NULL '
                    'THEN "' + null_syb + '" ELSE b.Expressed END, '
                    '"' + sep + '") AS Label '
                    'FROM '
                    '(SELECT a.Gene, b.Type FROM ' + self.name + '_ids AS a, '
                    ' ' + self.name + '_tissues AS b) AS a '
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
