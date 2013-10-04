'''
Handles the SQL connection for PAPPI.
The SQL is implemented here to use SQLite.
Everything else uses a standard SQLConnection interface,
such that the underlying database can easily be exchanged here.

@author: Patrick Flick
'''

from config import PAPPI_SQLITE_DEFAULT_DB
import sqlite3
import csv

PAPPI_SQL_CONN = None


def execute_script(script_filename, sql_conn=PAPPI_SQL_CONN):
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


def import_csv(csv_filename, table, csv_delimiter, has_header,
               column_names=None, column_types=None,
               csv_quoting=None, sql_conn=PAPPI_SQL_CONN):
    """
    Imports a CSV file into an SQL table.

    Parameters:
        - csv_filename:     The file name of the CSV file to be imported.
        - table:            The name of the SQL table to be created.
        - csv_delimiter:    The delimiter/seperator of the CSV file.
        - has_header:       Whether the CSV file has a header row.
        - column_names:     Names for the columns, if none are given either the
                            header row of the CSV file is used for column names
                            or general names are given: Column_i with i={1,..}
        - column_types:     The SQL types of the columns. Default: varchar(16)
        - csv_quoting:      Field-Quoting of the CSV file. Default: QUOTE_NONE
        - sql_conn:         The SQL connection to be used. Default: current
                            connection.
    """
    # check parameters
    if not csv_quoting:
        csv_quoting = csv.QUOTE_NONE

    # open the input file
    with open(csv_filename, 'r') as csv_file:
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
                row = csv_reader.next()
                for item in row:
                    # replace spaces with underscores and add to column names
                    column_names.append(item.replace(" ", "_"))
            else:
                row = csv_reader.next()
                for i in range(1, len(row) + 1):
                    column_names.append("Column_" + str(i))
                # reset the csv_reader object
                csv_reader = csv.reader(csv_file, delimiter=csv_delimiter,
                                        quoting=csv_quoting)
        else:
            if has_header:
                # ignore header
                csv_reader.next()

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
        cur.executemany('INSERT INTO "' + table + '" VALUES (' + vals + ')',
                        csv_reader)

        # close cursor and commit
        cur.close()
        sql_conn.commit()


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


def get_conn(db=PAPPI_SQLITE_DEFAULT_DB):
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
