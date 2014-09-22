'''
Handles the SQL connection for PAPPI.
The SQL is implemented here to use SQLite.
Everything else uses a standard SQLConnection interface,
such that the underlying database can easily be exchanged here.

@author: Patrick Flick
'''

from .data_config import DATABASE
import sqlite3
import csv
import re

PAPPI_SQL_CONN = None


def execute_script(script_filename, sql_conn=PAPPI_SQL_CONN):
    """
    Executes the given file as a script for the given SQL connection.

    @param script_filename: The file name of the script file to be executed.
    @param sql_conn:        The SQL connection to be used.
    """
    with open(script_filename, 'r') as script_file:
        # initialize the cursor object
        cur = sql_conn.cursor()

        # read script
        sql_script = script_file.read()

        # execute the script
        cur.executescript(sql_script)

        # close cursor and commit
        cur.close()
        sql_conn.commit()


def new_table_from_query(new_table, query, sql_conn=PAPPI_SQL_CONN,
                         overwrite=True):
    """
    Executes the given query and creates a new table with name `new_table` from
    the output query. In case the table with the given name exists, the
    `overwrite` variable decides what happens. If it is set to `True`, the
    existing table is deleted/overwritten. Otherwise the table is not
    overwritten and nothing happens.

    @param new_table:       The name for the new table to be created.
    @param query:           The SQL query to create the new table with.
    @param sql_conn:        The SQL connection to be used for the operation.
    @param overwrite:       Whether or not an existing table with name
                            `new_table` is to be overwritten.
    """
    # NOTE: this function is just for conveniance, because the pattern
    # implemented here existed all over the project. This way only one
    # convenient call has to be made.

    # get SQl cursor
    cur = sql_conn.cursor()

    # depending on the `overwrite` parameter: drop table and set correct
    # CREATE statement
    if overwrite:
        cur.execute('DROP TABLE IF EXISTS ' + new_table)
        sql_prefix = 'CREATE TABLE "' + new_table + '" AS '
    else:
        sql_prefix = 'CREATE TABLE IF NOT EXISTS "' + new_table + '" AS '

    # execute the query and save in the new table
    cur.execute(sql_prefix + query)

    # close cursor and commit to server
    cur.close()
    sql_conn.commit()


def extend_row_iterator(base_iterator, num, indices=None):
    for row in base_iterator:
        if not indices is None:
            row = [row[i] for i in indices]
        if len(row) == num:
            yield row
        else:
            row += [""]*(num - len(row))
            yield row
    return


def import_csv(csv_filename, table, csv_delimiter, has_header,
               import_columns=None,
               column_names=None, column_types=None,
               csv_quoting=None, sql_conn=PAPPI_SQL_CONN, skip_rows=0,
               row_iterator_wrapper=None):
    """
    Imports a CSV file into an SQL table.

    @param csv_filename:    The file name of the CSV file to be imported.
    @param table:           The name of the SQL table to be created.
    @param csv_delimiter:   The delimiter/seperator of the CSV file.
    @param has_header:      Whether the CSV file has a header row.
    @param import_columns:  A list of column indeces of columns to import into
                            the SQL database. If this is not set (=None), then
                            by default ALL columns are imported.
    @param column_names:    Names for the columns, if none are given either the
                            header row of the CSV file is used for column names
                            or general names are given: Column_i with i={1,..}
    @param column_types:    The SQL types of the columns. Default: varchar(16)
    @param csv_quoting:     Field-Quoting of the CSV file. Default: QUOTE_NONE
    @param sql_conn:        The SQL connection to be used. Default: current
                            connection.
    """
    # check parameters
    if not csv_quoting:
        csv_quoting = csv.QUOTE_NONE

    # open the input file
    with open(csv_filename, 'r') as csv_file:
        # in case rows are supposed to be skipped -> skip them:
        if skip_rows > 0:
            for i in range(0, skip_rows):
                l = csv_file.readline()

        # initialize the cursor object
        cur = sql_conn.cursor()

        # get csv reader for the CCSB file
        csv_reader = csv.reader(csv_file, delimiter=csv_delimiter,
                                quoting=csv_quoting)

        # if no column names are given, use the header row if available
        # or use default names Column_1, Column_2, ...
        if not column_names:
            column_names = []
            if has_header:
                # read header
                row = csv_reader.__next__()
                for item in row:
                    # replace spaces with underscores and add to column names
                    col_name = re.sub('[^0-9a-zA-Z_]+', '_', item)
                    column_names.append(col_name)
            else:
                row = csv_reader.__next__()
                for i in range(1, len(row) + 1):
                    column_names.append("Column_" + str(i))
                # reset the csv_reader object
                csv_reader = csv.reader(csv_file, delimiter=csv_delimiter,
                                        quoting=csv_quoting)
            # use only those wanted
            if not import_columns is None:
                column_names = [column_names[i] for i in import_columns]
        else:
            if has_header:
                # ignore header
                csv_reader.__next__()

        # if no column types are given, use the default for all
        default_col_type = "varchar(16)"
        if not column_types:
            column_types = [default_col_type] * len(column_names)

        # get sql column name and type string:
        cols = zip(column_names, column_types)
        sql_table_def = ", ".join('"' + x + '" ' + y for x, y in cols)

        # delete old and create new table
        cur.execute('DROP TABLE IF EXISTS "' + table + '"')
        cur.execute('CREATE TABLE "' + table + '" (' + sql_table_def + ')')

        # insert all lines
        vals = ", ".join(['?'] * len(column_names))
        if row_iterator_wrapper is None:
            row_iter = csv_reader
        else:
            row_iter = row_iterator_wrapper(csv_reader)
        cur.executemany('INSERT INTO "' + table + '" VALUES (' + vals + ')',
                        extend_row_iterator(row_iter, len(column_names),
                                            import_columns))

        # close cursor and commit
        cur.close()
        sql_conn.commit()


def linearize_table(src_table, excl_columns, cat_col_name, val_col_name,
                    dst_table, con=PAPPI_SQL_CONN):
    col_names = get_column_names(src_table, con)

    transpose_cols = [x for x in col_names if not x in excl_columns]

    remaining_cols = ", ".join(excl_columns)

    select_stmts = ['SELECT ' + remaining_cols + ' , \'' + x + '\' AS '
                    + cat_col_name + ', [' + x + '] AS ' + val_col_name
                    + ' FROM ' + src_table for x in transpose_cols]
    union = " UNION ".join(select_stmts)

    # get SQL cursor
    cur = con.cursor()
    cur.execute('DROP TABLE IF EXISTS ' + dst_table)
    cur.execute('CREATE TABLE ' + dst_table + ' AS ' + union)

    # commit SQl changes
    cur.close()
    con.commit()


def dump_csv(outfile, table, sql_conn=PAPPI_SQL_CONN):
    """
    Dumps the given SQL table as CSV file into the given output file.

    @param outfile:  The file to be written to.
    @param table:    The SQL table to be dumped.
    @param sql_conn: The SQL connection to be used.
    """
    # get SQL curser object
    cur = sql_conn.cursor()

    cols = []
    for row in cur.execute('PRAGMA table_info("' + table + '")'):
        cols.append(row[1])

    wr = csv.writer(outfile, delimiter=' ', quoting=csv.QUOTE_NONE)
    wr.writerow(cols)

    for row in cur.execute('SELECT * FROM "' + table + '"'):
        wr.writerow(row)

    cur.close()
    sql_conn.commit()


def table_exists(table, sql_conn=PAPPI_SQL_CONN):
    """
    Returns whether a table with the given name exists in the SQL database
    given by the SQL connection.

    @param table:       The name of the SQL table to be checked for.
    @param sql_conn:    The SQl Connection object to be used.
                        Default: current connection.
    @return             Boolean, whether the table with the given name exits.
    """
    cur = sql_conn.cursor()
    # the general SQL version
    #cur.execute('SELECT COUNT(*) FROM INFORMATION_SCHEMA.TABLES '
                #'WHERE TABLE_NAME = "' + table + '"')
    # for SQLite:
    cur.execute('SELECT COUNT(*) FROM sqlite_master '
                'WHERE type="table" AND name="' + table + '"')
    result = cur.fetchone()[0]
    return result > 0


def get_column_names(table, sql_conn=PAPPI_SQL_CONN):
    """
    Returns the names of all columns in the given table.

    @param table:       The table which columns are to be returned.
    @param sql_conn:    The SQL connection object.
    @Returns            A list of names of the columns.
    """
    # source for this:
    # http://stackoverflow.com/questions/7831371/is-there-a-way-to-get-a-list-of-column-names-in-sqlite
    cur = sql_conn.cursor()
    cur.execute('SELECT * FROM ' + table)
    names = list(map(lambda x: x[0], cur.description))
    cur.close()
    # return the result
    return names


def get_conn(db=DATABASE):
    """
    Returns the current SQL connection object. Connects and returns a
    new object if there is no current SQL connection.

    @param db: The SQLite database name to be used (default: pappiDB.sqlite)
    @return: The current SQL connection.
    """
    global PAPPI_SQL_CONN
    if (not PAPPI_SQL_CONN):
        # establish a SQLite3 connection.
        # If another database is to be used, this is the only
        # place where changes are necessary
        PAPPI_SQL_CONN = sqlite3.Connection(db)
    return PAPPI_SQL_CONN


def close(con):
    """
    Closes the active connection.
    """
    global PAPPI_SQL_CONN
    if (PAPPI_SQL_CONN):
        PAPPI_SQL_CONN.close()
        PAPPI_SQL_CONN = None
