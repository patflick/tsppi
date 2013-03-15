import sqlite3
import csv
import os

from . import PAPPI_SQL_FOLDER_PATH

PAPPI_DEFAULT_SQLITE_DB = 'pappiDB.sqlite'
PAPPI_HPA_RAW_TABLE_NAME = 'hpa_normal_tissue'
PAPPI_HPA_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'hpa_filter.sql')

def init_hpa_data(database):
    """
    Initializes the imported HPA data. This includes filtering for
    APE scoring with `High` or `Medium` Reliability and adding
    a numerical value for the Level (High: 3, Medium: 2, Low: 1, None: 0).
    """
    with open(PAPPI_HPA_FILTER_SCRIPT, 'r') as script_file:
        # initialize the sqlite3 connection
        con = sqlite3.Connection(database)
        cur = con.cursor()

        # read script
        sql_script = script_file.read()
        
        # execute the script
        cur.executescript(sql_script)

        # commit and close sql connection
        cur.close()
        con.commit()
        con.close()


def import_hpa_file(infile, database, table = PAPPI_HPA_RAW_TABLE_NAME):
    # initialize the sqlite3 connection
    con = sqlite3.Connection(database)
    cur = con.cursor()
    
    # create table for the raw HPA data:
    # "Gene","Tissue","Cell type","Level","Expression type","Reliability"
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("Gene" varchar(16), "Tissue" varchar(256), "Cell.type" varchar(256), "Level" varchar(32), "Expression.type" varchar(12), "Reliability" varchar(32));')
    
    # get csv reader for the hpa file
    csv_reader = csv.reader(infile, delimiter=',',quoting=csv.QUOTE_ALL)
    # ignore header line
    csv_reader.next()
    
    # insert all lines
    cur.executemany('INSERT INTO "' + table + '" VALUES (?, ?, ?, ?, ?, ?)', csv_reader)
    
    # commit and close sql connection
    cur.close()
    con.commit()
    con.close()
    
def import_hpa(hpafile, database=PAPPI_DEFAULT_SQLITE_DB):
    # first import the file as table
    import_hpa_file(hpafile, database, PAPPI_HPA_RAW_TABLE_NAME)
    # then create/fill the other tables (scoring, mean tissue, summary)
    init_hpa_data(database)