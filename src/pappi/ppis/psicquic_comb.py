from . import ppi
from .. import sql
import re


class PsicquicAll(ppi.PPI):
    """
    Creates a combined PPI network from all PSICQUIC networks
    """

    def __init__(self, sql_connection, ppi_names):
        """
        Creates a new PPI network by combining the networks from `ppi_names`.

        Parameter:
            - sql_connection:   SQL connection object for the database.
            - ppi_names:        Names of the PSICQUIC PPI networks in the
                                datebase to be combined.
        """
        # copy references into the new instance
        ppi.PPI.__init__(self, "", sql_connection, "psicquic_all",
                         "Gene1", "Gene2", "hgnc")
        self.src_ppis = ppi_names

    def init_ppi(self, verbose=False):
        """
        Creates the combined PPI network as a new database table.
        """
        union = " UNION ALL ".join("SELECT Gene1, Gene2 FROM psicquic_" + name
                                   for name in self.src_ppis)
        sqlquery = "SELECT * FROM (" + union + ")"
        sql.new_table_from_query(self.name + "_union", sqlquery, self.sql_conn)
        sqlquery = ("SELECT Gene1, Gene2, COUNT(*) AS Count FROM " + self.name
                    + "_union GROUP BY Gene1, Gene2")
        sql.new_table_from_query(self.name + "_count", sqlquery, self.sql_conn)
        sqlquery = ("SELECT Gene1, Gene2 FROM " + self.name + "_count")
        sql.new_table_from_query(self.name, sqlquery, self.sql_conn)
