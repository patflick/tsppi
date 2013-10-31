class TableManager:
    """
    A class to handle import of csv files into SQL tables and
    handle multiple sub-steps of processing these SQL tables.
    This class is a super-class to both the ppis and the
    expression data sets.
    """
    tmp_table_idx = 0
    cur_tmp_name = None

    def __init__(self, table_name, sql_connection):
        self.name = table_name
        self.sql_conn = sql_connection

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
        if self.cur_tmp_name:
            return self.cur_tmp_name
        else:
            return next_tmp_table()

    def next_tmp_table(self, suffix=None):
        """
        Creates a new temporary SQL table name and returns that new name.
        The next call to get_cur_tmp_table() will return this name.

        All subclasses MUST NOT overwrite this function and must make
        use of it to get the names for temporary tables.
        """
        self.tmp_table_idx += 1
        if suffix is None:
            suffix = str(self.tmp_table_idx)

        # set new table name
        if suffix == "":
            self.cur_tmp_name = self.name
        else:
            self.cur_tmp_name = self.name + '_' + suffix

        # return via the current table function
        return self.get_cur_tmp_table()
