from .expr import GeneExpression
from .. import id_mapping
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
        # TODO table naming scheme as in ppis (maybe create another
        #       super-super class to contain the table management!?)
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table()
        id_mapping.map_identifier(src_table, ["Gene"], "ensembl", dst_table,
                                  id_mapping.UNIFYING_ID, self.sql_conn, True)
