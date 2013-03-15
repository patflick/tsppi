
import sqlite3
import csv
import re
import os

from . import matching
from . import PAPPI_SQL_FOLDER_PATH

# TODO put all constants in an own module
PAPPI_PPI_RAW_TABLE_NAME = 'ppi_raw_proteins'
PAPPI_PPI_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'ppi_filter.sql')
PAPPI_DEFAULT_SQLITE_DB = 'pappiDB.sqlite'
PAPPI_ENSP2ENSG_TABLE_NAME = 'ensg_to_ensp'

def import_stringdb_file(infile, database, table=PAPPI_PPI_RAW_TABLE_NAME):
    # initialize the sqlite3 connection
    con = sqlite3.Connection(database)
    cur = con.cursor()
    
    # create table for the raw HPA data:
    # "Gene","Tissue","Cell type","Level","Expression type","Reliability"
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, "Protein1" varchar(16), "Protein2" varchar(16), "Score" int)')
    
    # get csv reader for the hpa file
    csv_reader = csv.reader(infile, delimiter=' ',quoting=csv.QUOTE_NONE)
    
    # compile regex for extracting proteins (don't put the 9606. in the db)
    p = re.compile(r'^9606\.(ENSP\d+)$')

    # insert all lines
    for row in csv_reader:
        # extract the ensembl protein ID
        p1 = p.match(row[0]).group(1)
        p2 = p.match(row[1]).group(1)
        # insert into table
        cur.execute('INSERT INTO "' + table + '" (Protein1, Protein2, Score) VALUES (?, ?, ?)', (p1,p2,row[2]))

    # commit and close sql connection
    cur.close()
    con.commit()
    con.close()

def init_stringdb_ppi(database):
    with open(PAPPI_PPI_FILTER_SCRIPT, 'r') as script_file:
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

def import_stringdb(ppifile, ensg2ensp_file, database=PAPPI_DEFAULT_SQLITE_DB):
    # first import the file as table
    import_stringdb_file(ppifile, database, PAPPI_PPI_RAW_TABLE_NAME)
    
    # import ENSG<->ENSP matching
    matching.import_p2g_file(ensg2ensp_file, database, PAPPI_ENSP2ENSG_TABLE_NAME)
    
    # then create/fill the other tables (scoring, mean tissue, summary)
    init_stringdb_ppi(database)
