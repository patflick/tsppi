import csv
from .. import sql


class PPI:
    """ An interface to all PPI networks. """

    # default parameters for the raw file
    file_has_header = True
    file_field_seperator = '\t'
    file_quoting = csv.QUOTE_NONE
    ppi_name = "no_name"

    def __init__(self, filename, sql_connection):
        """
        Creates a new PPI network, with the given filename and sql
        connection object.

        Parameters:
            - filename:         Filename of the ppi network file.
            - sql_connection:   SQL connection object for the database,
                                into which the PPI network is to be
                                imported.
        """
        # copy references into the new instance
        self.filename = filename
        self.sql_conn = sql_connection

    def import_raw_file(self):
        """
        Imports a PPI file in csv format into the SQL database.
        """
        table_name = self.ppi_name + "_raw"
        sql.import_csv(self.filename, table_name, self.file_field_seperator,
                       self.file_has_header, csv_quoting=self.file_quoting)

    def filter(self):
        """
        Filters the PPI network for interactions with a certain confidence
        level. This is only implemented for some PPIs. E.g. string-db and
        PSICQUIC have these confidence scores by which the interactions can
        be filtered.
        """
        # the interface/base-class does not implement filtering.
        pass

    def id_mapping(self):
        """
        Maps the gene/protein identifiers to a unifying identfier.
        """
        # TODO is it possible to call a general mapping function here twice?
