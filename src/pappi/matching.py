import sqlite3
import csv

PAPPI_ENSP2ENSG_TABLE_NAME = 'ensg_to_ensp'

def import_p2g_file(infile, database, table=PAPPI_ENSP2ENSG_TABLE_NAME):
    # initialize the sqlite3 connection
    con = sqlite3.Connection(database)
    cur = con.cursor()
    
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
    
    # commit and close sql connection
    cur.close()
    con.commit()
    con.close()