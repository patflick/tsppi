'''
Import operations for ID matching tables. (Ensembl ENSG<->ENSP and Ensembl ENSG <-> EntreZ Gene)

@author: Patrick Flick
'''

import csv

PAPPI_ENSP2ENSG_TABLE_NAME = 'ensg_to_ensp'

def import_p2g_file(infile, sql_conn, table=PAPPI_ENSP2ENSG_TABLE_NAME):
    """
    Imports the Ensembl protein to gene matching ENSG<->ENSP.
    
    Get the ensembl ID matching from
    http://www.ensembl.org/biomart/martview/
    using homo sapiens version 70 database and output the attributes:
        - Ensembl Gene ID
        - Ensembl Protein ID
    with NO filters
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the data is to be imported.
    """
    # initialize the cursor object
    cur = sql_conn.cursor()
    
    # create table for the raw HPA data:
    # "Gene","Tissue","Cell type","Level","Expression type","Reliability"
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("ENSG" varchar(16), "ENSP" varchar(256));')
    
    # get csv reader for the hpa file
    csv_reader = csv.reader(infile, delimiter=',',quoting=csv.QUOTE_ALL)
    # ignore header line
    csv_reader.next()
    
    # insert all lines
    cur.executemany('INSERT INTO "' + table + '" VALUES (?, ?)', csv_reader)
    
    # close cursor and commit
    cur.close()
    sql_conn.commit()