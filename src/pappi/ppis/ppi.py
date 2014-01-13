import csv
from .. import sql
from .. import id_mapping
from ..table_manager import TableManager

# import the UNIFYING ID variable
from ..id_mapping import UNIFYING_ID


class PPI(TableManager):
    """ An interface to all PPI networks. """

    # default parameters for the raw file
    file_has_header = True
    file_field_seperator = '\t'
    file_quoting = csv.QUOTE_NONE

    def __init__(self, filename, sql_connection, name,
                 gene1_colname, gene2_colname, orig_id):
        """
        Creates a new PPI network. This MUST explicitly be called by all
        subclasses of this class to properly initialize all class members.

        Parameters:
            - filename:         Filename of the ppi network file.
            - sql_connection:   SQL connection object for the database,
                                into which the PPI network is to be
                                imported.
            - name:             The name for the PPI network. This is
                                used for example for the names of the
                                SQL tables.
            - gene1_colname:    The name of the column of the first gene
                                of the interaction.
            - gene2_colname:    The name of the column of the second gene
                                of the interaction.
            - orig_id:          The gene identifier type (e.g. 'uniprot',
                                'ensembl', ...) of the gene columns.
        """
        # call constructor of super-class
        TableManager.__init__(self, name, sql_connection)
        # copy references into the new instance
        self.filename = filename
        self.gene1_colname = gene1_colname
        self.gene2_colname = gene2_colname
        self.orig_id = orig_id

    def init_ppi(self, verbose=False):
        """
        Performs all the steps of initialization of the PPI network
        including loading the raw file into the SQl database,
        filtering the interactions, mapping to a unifying id and
        normalizing the PPI network to a common format.
        """
        if verbose:
            print("Initializing PPI '" + self.name + "'")
            print("    importing raw file")
        self.import_raw_file()
        if verbose:
            print("    filtering PPI")
        self.filter()
        if verbose:
            print("    mapping ids")
        self.id_mapping()
        if verbose:
            print("    normalizing table")
        self.normalize_table()
        if verbose:
            print("    normalizing graph")
        self.normalize_graph()

    def import_raw_file(self):
        """
        Imports a PPI file in csv format into the SQL database.
        """
        table_name = self.next_tmp_table("raw")
        sql.import_csv(self.filename, table_name, self.file_field_seperator,
                       self.file_has_header, csv_quoting=self.file_quoting,
                       sql_conn=self.sql_conn)

    def filter(self):
        """
        Filters the PPI network for interactions with a certain confidence
        level. This is only implemented for some PPIs. E.g. string-db and
        PSICQUIC have these confidence scores by which the interactions can
        be filtered.
        """
        # the interface/base-class does not implement filtering.
        pass

    def id_mapping(self):
        """
        Maps the gene/protein identifiers to a unifying identfier.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table()
        id_mapping.map_identifier(src_table, [self.gene1_colname,
                                  self.gene2_colname], self.orig_id,
                                  dst_table, UNIFYING_ID, self.sql_conn,
                                  verbose=True)

    def normalize_table(self):
        """
        Maps the columns of the original PPI table
        to columns of names: Gene1, Gene2. So that all tables have
        the same format.
        """
        # this happens right after the id mapping, just take the two mapped
        # ID columns and put them into a new table with column names
        # Gene1 and Gene2
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('normalized')
        sqlquery = ('SELECT '
                    '' + self.gene1_colname + ' AS Gene1, '
                    '' + self.gene2_colname + ' AS Gene2 '
                    'FROM ' + src_table)
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def normalize_graph(self):
        """
        Normalizes the PPI network by making it undirected with one
        row (i.e. pair of proteins) per edge (no back and forth edges).
        """
        # this works on the normalized columns with
        # the same identifier for all PPIs, thus
        # this method should not be overwritten
        # by subclasses
        src_table = self.get_cur_tmp_table()
        # this is the last step in the PPI pre-processing,
        # thus the result does not get put into a temporary table
        dst_table = self.next_tmp_table("")
        sqlquery = ('SELECT DISTINCT '
                    'CASE WHEN Gene1 > Gene2 '
                    '  THEN Gene2 ELSE Gene1 '
                    'END AS Gene1, '
                    'CASE WHEN Gene1 > Gene2 '
                    '  THEN Gene1 ELSE Gene2 '
                    'END AS Gene2 '
                    'FROM ' + src_table)
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def create_all_ids_table(self):
        """
        Creates a table named `ppi_name`_ids that holds all the distinct IDs
        used in the ppi network.
        """
        cur = self.sql_conn.cursor()
        if not sql.table_exists(self.name + '_ids', self.sql_conn):
            cur.execute('CREATE TABLE IF NOT EXISTS `' + self.name + '_ids` '
                        '(id integer primary key autoincrement, '
                        'Gene varchar(16))')
            sqlquery = ('SELECT Gene1 as Gene FROM ' + self.name + ' UNION '
                        'SELECT Gene2 as Gene FROM ' + self.name)
            cur.execute('INSERT INTO `' + self.name + '_ids` (Gene) '
                        + sqlquery)
        # create an index for fast merging
        cur.execute('DROP INDEX IF EXISTS ' + self.name + '_ids_index')
        cur.execute('CREATE UNIQUE INDEX ' + self.name + '_ids_index '
                    'ON ' + self.name + '_ids (Gene)')
        cur.close()
        self.sql_conn.commit()

    def export_to_edge_list(self, only_ids=None):
        """
        Exports the network from the SQL table into an edge list of IDs
        in the range [1,#IDs] and saves the resulting list into a tab
        separated file.

        @param filename:        The file-path for the file to which the
                                edgelist is written.
        @param only_ids:        An SQL table of Gene ids. The exported network
                                will consist only of edges between those ids,
                                i.e. the subnetwork/subgraph defined by this
                                set of IDs. If this is set to `None`, all
                                edges are exported. (default: None)
        """
        # check if ids table exists (if not -> create it)
        if not sql.table_exists(self.name + '_ids', self.sql_conn):
            self.create_all_ids_table()

        # set the id default table to be used for mapping the Gene names to
        # continuous integers
        id_table = self.name + '_ids'

        # check if the ID table has to be filtered
        if not only_ids is None:
            sqlquery = ('SELECT * FROM ' + id_table + ' WHERE Gene IN '
                        '(SELECT DISTINCT Gene FROM ' + only_ids + ')')
            # either give as inner query or create table and return new table
            # name
            id_table = '(' + sqlquery + ')'

        # merge PPI network  with unique integer ids (node-ids)
        sqlquery = ('SELECT b.id AS Gene1, c.id AS Gene2 '
                    'FROM ' + self.name + ' AS a '
                    'INNER JOIN ' + id_table + ' AS b '
                    ' ON a.Gene1 = b.Gene '
                    'INNER JOIN ' + id_table + ' AS c '
                    ' ON a.Gene2 = c.Gene ')

        sql.new_table_from_query(self.name + "_for_export", sqlquery,
                                 self.sql_conn)

        # save results to file
        #with open(filename, "w") as f:
        #    for row in cur.fetchall():
        #        f.write("%i\t%i\n" % (row[0], row[1]))
