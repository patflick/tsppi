from .expr import GeneExpression
from .. import id_mapping
from .. import sql
import csv


class Emtab(GeneExpression):

    def __init__(self, filename, sql_connection):
        """
        """
        # firstly initialize the superclass
        GeneExpression.__init__(self, filename, sql_connection, "emtab")

        # set the raw file properties
        self.file_field_seperator = "\t"

    def import_raw_file(self):
        """
        Imports the expression file in raw format into the SQL database

        For the E-MTAB-513 expression file, first the 3 commented lines
        have to be removed. Thus this function overwrites the base
        class function.
        """
        table_name = self.next_tmp_table("raw")
        sql.import_csv(self.filename, table_name, self.file_field_seperator,
                       self.file_has_header, csv_quoting=self.file_quoting,
                       sql_conn=self.sql_conn, skip_rows=3)

    def linearize_table(self):
        """
        The E-MTAB-513 table is in a [Gene]x[Tissue] format, this has to
        be linearized before continuing.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('linear')
        sql.linearize_table(src_table, ["Gene_ID", "Gene_Name"], "Tissue",
                            "ExpressionValue", dst_table, self.sql_conn)

    def id_mapping(self):
        """
        Maps the gene identifier to the unifying ID defined in `id_mapping`.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table()
        id_mapping.map_identifier(src_table, ["Gene_Name"], "hgnc", dst_table,
                                  id_mapping.UNIFYING_ID, self.sql_conn, True)

    # E-MTAB-513 doesn't supply reliablilites, thus no filter operation will
    # be performed

    def normalize_table(self):
        """
        Puts the table into the format: [Gene], [Tissue], [ExpressionValue]
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('normalized')
        sqlquery = ('SELECT Gene_Name as Gene, Tissue AS Type, '
                    'Tissue, ExpressionValue '
                    'FROM ' + src_table)
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)

    def classify(self):
        """
        Classifies the ExpressionValue into a binary value Expressed.
        """
        # set the classification condition
        self.classify_cond = ">= 1"
        # call the super-class function
        GeneExpression.classify(self)
