from .expr import GeneExpression
from .hpa import HPA
from .. import id_mapping
from .. import sql
import csv


class HPA_All(HPA):
    """
    A sub class from HPA which does not filter out any data based on
    reliability (in contrast to the regular HPA class, which only
    keeps reliable data in the pipeline)
    """

    def __init__(self, filename, sql_connection):
        """
        Constructs the Human Protein Atlas gene expression dataset
        """
        # firstly initialize the HPA superclass
        HPA.__init__(self, filename, sql_connection)
        # set the table name differently
        self.name = "hpa_all"

        # set the raw file properties
        self.file_quoting = csv.QUOTE_ALL
        self.file_field_seperator = ","

    def filter(self):
        """
        In the HPA_All class, nothing is filtered out (in contrast to the
        regular HPA class, which only keeps reliable data).
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table('filtered')
        sqlquery = ('SELECT Gene, Tissue, Cell_Type, Level, '
                    'Expression_type '
                    'FROM ' + src_table)
        sql.new_table_from_query(dst_table, sqlquery, self.sql_conn)
