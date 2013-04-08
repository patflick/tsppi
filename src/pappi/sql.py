'''
Handles the SQL connection for PAPPI.
The SQL is implemented here to use SQLite.
Everything else uses a standard SQLConnection interface,
such that the underlying database can easily be exchanged here.

@author: Patrick Flick
'''

from config import PAPPI_SQLITE_DEFAULT_DB
import sqlite3

PAPPI_SQL_CONN=None


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
    if (PAPPI_SQL_CONN):
        PAPPI_SQL_CONN.close()
        PAPPI_SQL_CONN=None