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

    def export_node_labels(self, node_ids_tbl):
        """
        Exports binary expression node labels for the given Node-IDs table
        which must be of the form (id, Gene).
        """
        self.create_tissue_table()
