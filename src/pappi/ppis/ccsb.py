from . import ppi


class CCSB(ppi.PPI):
    """
    The PPI network HI-2012-Pre from CCSB at
    Dana Farber / Harvard Medical School
    http://interactome.dfci.harvard.edu/H_sapiens/index.php?page=download
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
        # copy references into the new instance
        ppi.PPI.__init__(self, filename, sql_connection, "ccsb",
                         "GENE_IDA", "GENE_IDB", "entrez")

        # init custom fields
        self.has_header = True
        self.file_field_seperator = '\t'
