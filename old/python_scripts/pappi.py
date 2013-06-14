import sqlite3
import csv
import re

# the script must be in the same folder
PAPPI_HPA_FILTER_SCRIPT = 'hpa_filter.sql'
PAPPI_PPI_FILTER_SCRIPT = 'ppi_filter.sql'

PAPPI_DEFAULT_SQLITE_DB = 'pappiDB.sqlite'
PAPPI_HPA_RAW_TABLE_NAME = 'hpa_normal_tissue'
PAPPI_PPI_RAW_TABLE_NAME = 'ppi_raw_proteins'
PAPPI_ENSP2ENSG_TABLE_NAME = 'ensg_to_ensp'


def init_hpa_data(database):
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


def import_ppi_file(infile, database, table=PAPPI_PPI_RAW_TABLE_NAME):
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

def init_ppi(database):
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

def import_ppi(ppifile, ensg2ensp_file, database=PAPPI_DEFAULT_SQLITE_DB):
    # first import the file as table
    import_ppi_file(ppifile, database, PAPPI_PPI_RAW_TABLE_NAME)
    
    # import ENSG<->ENSP matching
    import_ppi_file(ensg2ensp_file, database, PAPPI_ENSP2ENSG_TABLE_NAME)
    
    # then create/fill the other tables (scoring, mean tissue, summary)
    init_ppi(database)



