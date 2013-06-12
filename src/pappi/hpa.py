'''
Import and filtering operations for the HPA data.

@author: Patrick Flick
'''

import csv

from config import PAPPI_SQL_HPA_FILTER_SCRIPT
from config import PAPPI_SQL_HPA_GENE_LEVELS_SCRIPT
from sql import execute_script


PAPPI_HPA_RAW_TABLE_NAME = 'hpa_normal_tissue'


def init_tissue_data(sql_conn):
    """
    Initializes the imported HPA data. This includes filtering for
    APE scoring with `High` or `Medium` Reliability and adding
    a numerical value for the Level (High: 3, Medium: 2, Low: 1, None: 0).
    
    @param sql_conn: The SQL connection to be used.
    """
    execute_script(PAPPI_SQL_HPA_FILTER_SCRIPT, sql_conn)


def init_gene_levels(sql_conn):
    """
    Creates a new `hpa_gene_levels` table with summarizing
    counts for all expression levels per gene.
    
    @param sql_conn: The SQL connection to be used.
    """
    execute_script(PAPPI_SQL_HPA_GENE_LEVELS_SCRIPT, sql_conn)



def import_tissue_file(infile, sql_conn, table = PAPPI_HPA_RAW_TABLE_NAME):
    """
    Imports the downloaded HPA tissue file into the database given by the
    `sql_conn` SQL connection into a new table given by `table`.
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the HPA data is to be imported.
    """
    # initialize the cursor object
    cur = sql_conn.cursor()
    
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
    
    # close cursor and commit
    cur.close()
    sql_conn.commit()



def import_tissue(hpafile, sql_conn):
    """
    Imports the HPA data, filters and initializes further tables (by calling init_tissue_data).
    
    @param hpafile: The file handle for the HPA tissue data file.
    @param sql_conn: The SQL connection to be used.
    """
    # first import the file as table
    import_tissue_file(hpafile, sql_conn, PAPPI_HPA_RAW_TABLE_NAME)
    # then create/fill the other tables (scoring, mean tissue, summary)
    init_tissue_data(sql_conn)
