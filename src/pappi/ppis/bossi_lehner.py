# import the module for the super-class PPI
from .ppi import PPI


class Bossi_Lehner(PPI):
    """
    The PPI network from the Bossi and Lehner paper:
    http://www.ncbi.nlm.nih.gov/pubmed/19357639

    This is a combined network from multiple sources.
    """

    def __init__(self, filename, sql_connection):
        """
        Creates a new PPI network, with the given filename and sql
        connection object.

        Parameter:
            - filename:         Filename of the ppi network file.
            - sql_connection:   SQL connection object for the database,
                                into which the PPI network is to be
                                imported.
        """
        # create super class with all necessary paramenters
        PPI.__init__(self, filename, sql_connection, "bossi",
                     "Gene1", "Gene2", "ensembl")

        self.has_header = True
        self.file_field_seperator = '\t'

    def filter(self):
        """
        Overwrites the filter method in order to throw
        out all the expression data from the raw data
        and just be left with the interaction data set.
        """
        src_table = self.get_cur_tmp_table()
        dst_table = self.next_tmp_table()
        cur = self.sql_conn.cursor()
        cur.execute('CREATE TABLE IF NOT EXISTS ' + dst_table + ' AS '
                    'SELECT Gene1, Gene2 FROM ' + src_table)

        cur.close()
        self.sql_conn.commit()
