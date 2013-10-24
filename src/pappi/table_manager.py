class TableManager:
    """
    A class to handle import of csv files into SQL tables and
    handle multiple sub-steps of processing these SQL tables.
    This class is a super-class to both the ppis and the
    expression data sets.
    """
    tmp_table_idx = 0

    def __init__(self, table_name, sql_connection):
        self.name = table_name
        self.sql_conn = sql_connection
        # set the table index WHAT TAH FAAACK?
        self.tmp_table_idx = 42

    def get_cur_tmp_table(self):
        """
        Returns the name of the temporary SQL table with the results of
        the current step in the pre-processing.

        This function is needed, because subclasses of PPI might decide
        to skip a few pre-processing steps in case they are not necessary
        for that specific PPI. In that case fewer temporary tables
        are needed to complete the pre-processing.

        All subclasses MUST NOT overwrite this function and must make
        use of it to get the names for temporary tables.
        """
        if self.tmp_table_idx == 0:
            return self.name + "_raw"
        else:
            return self.name + "_tmp_" + str(self.tmp_table_idx)

    def next_tmp_table(self):
        """
        Creates a new temporary SQL table name and returns that new name.
        The next call to get_cur_tmp_table() will return this name.

        All subclasses MUST NOT overwrite this function and must make
        use of it to get the names for temporary tables.
        """
        self.tmp_table_idx += 1
        return self.get_cur_tmp_table()
