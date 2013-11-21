from . import ppi
from .. import sql


class StringDB(ppi.PPI):
    """
    Imports and filters the string-db PPI network.
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
        ppi.PPI.__init__(self, filename, sql_connection, "string",
                         "Gene1", "Gene2", "ensp")

        # init custom fields
        self.has_header = False
        self.file_field_seperator = ' '

    def import_raw_file(self):
        """
        Imports the string-db PPI file in csv format into the SQL database.
        """
        # custom row iterator wrapper:
        def string_row_iter(base_iter):
            for row in base_iter:
                row = [s.replace("9606.", "") for s in row]
                yield row
            return
        # import the string PPI without the 9606. prefix
        table_name = self.next_tmp_table("raw")
        sql.import_csv(self.filename, table_name, self.file_field_seperator,
                       self.file_has_header,
                       column_names=["Gene1", "Gene2", "Reliability"],
                       column_types=["varchar(16)", "varchar(16)", "int"],
                       row_iterator_wrapper=string_row_iter,
                       csv_quoting=self.file_quoting, sql_conn=self.sql_conn)

    def filter(self):
        # don't filter yet
        # TODO (maybe) implement other SW design to handle differently filtered
        #              PPIs
        pass
