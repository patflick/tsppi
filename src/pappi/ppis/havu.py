# import the module for the super-class PPI
from .ppi import PPI


class Havugimana(PPI):
    """
    The PPI network from the Havugimana paper:
    A census of human soluble protein complexes
    http://www.ncbi.nlm.nih.gov/pubmed/22939629
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
        PPI.__init__(self, filename, sql_connection, "havu",
                     "Gene1", "Gene2", "uniprot")

        self.has_header = True
        self.file_field_seperator = '\t'
