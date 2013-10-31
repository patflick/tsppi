from .expr import GeneExpression
from .. import id_mapping
from .. import sql
import csv


class HPA(GeneExpression):

    def __init__(self, filename, sql_connection):
        """
        Constructs the Human Protein Atlas gene expression dataset
        """
        # firstly initialize the superclass
        GeneExpression.__init__(self, filename, sql_connection, "hpa")

        # set the raw file properties
        self.file_quoting = csv.QUOTE_ALL
        self.file_field_seperator = ","

    def id_mapping(self):
        """
        Maps the gene identifier to the unifying ID defined in `id_mapping`.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table()
        id_mapping.map_identifier(src_table, ["Gene"], "ensembl", dst_table,
                                  id_mapping.UNIFYING_ID, self.sql_conn, True)

    def filter(self):
        """
        Filters out all unreliable values
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('filtered')
        sqlquery = ('SELECT Gene, Tissue, Cell_Type, Level, Expression_type '
                    'FROM ' + src_table + ' '
                    'WHERE (Expression_type = "APE" AND '
                    '(Reliability = "Medium" OR Reliability = "High")) '
                    'OR (Expression_type = "Staining" AND '
                    'Reliability = "Supportive")')
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def normalize_table(self):
        """
        Normalizes the table into the format [Gene], [Type],
        [ExpressionValue], ...
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('normalized')
        sqlquery = ('SELECT Gene, Tissue || " - " || Cell_Type AS Type, '
                    'Tissue, Cell_Type, Level AS ExpressionValue '
                    'FROM ' + src_table)
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def classify(self):
        """
        Classifies the ExpressionValue into expressed or non-expressed.
        """
        # set classification condition
        self.classify_cond = 'NOT IN ("None", "Negative")'
        # call super class function
        GeneExpression.classify(self)
