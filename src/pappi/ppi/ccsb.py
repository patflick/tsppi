import ppi


class CCSB(ppi.PPI):

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
        ppi.PPI.__init__(self, filename, sql_connection)

        # init custom fields
        self.ppi_name = "ccsb"
        self.has_header = True
        self.file_field_seperator = ' '

    def import_raw_file(self):
        # the super class implementation suffices, because the CSV file
        # format is very simple
        ppi.PPI.import_raw_file(self)

    def normalize(self):
        """
        Matches everything to HGNC ids and makes the PPI network
        undirected with one pair per edge (no back and forth edges).
        """

        # TODO: execute SQL statements here (no more outside SQL scripts)
        # TODO: 1st) ID Matching
        # TODO: 2nd) Normalization to one edge per pair!
