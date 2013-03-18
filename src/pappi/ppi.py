'''
Import and filtering operations for PPI data.

@author: Patrick Flick
'''

import csv
import re

from . import matching
from config import PAPPI_SQL_PPI_FILTER_SCRIPT

# TODO put all constants in an own module
PAPPI_PPI_RAW_TABLE_NAME = 'ppi_raw_proteins'

PAPPI_ENSP2ENSG_TABLE_NAME = 'ensg_to_ensp'


def import_stringdb_file(infile, sql_conn, table=PAPPI_PPI_RAW_TABLE_NAME):
    """
    Imports the string-db PPI network from the tab separated file. The PPI is imported
    using the `sql_conn` SQL connection into a new table given by `table`.
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the string-db data is to be imported.
    """
    # initialize the cursor object
    cur = sql_conn.cursor()
    
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

    # close cursor and commit
    cur.close()
    sql_conn.commit()



def init_stringdb_ppi(sql_conn):
    """
    Initializes the imported string-db PPI data. This includes replacing
    proteins (ENSP IDs) with genes (ENSG IDs) and filtering out duplicate edges.
    
    @param sql_conn: The SQL connection to be used.
    """
    with open(PAPPI_SQL_PPI_FILTER_SCRIPT, 'r') as script_file:
        # initialize the cursor object
        cur = sql_conn.cursor()

        # read script
        sql_script = script_file.read()
        
        # execute the script
        cur.executescript(sql_script)

        # close cursor and commit
        cur.close()
        sql_conn.commit()



def import_stringdb(ppifile, ensg2ensp_file, sql_conn):
    """
    Imports and initializes the string-db PPI network. Also the given ENSG<->ENSP mapping file
    is imported, which is used during the string-db initialization.
    
    @param ppifile: The file handle for the PPI network file, downloaded from string-db.org
    @param ensg2ensp_file: The file handle for the ENSG<->ENSP mapping to be imported.
    @param sql_conn: The SQL connection to be used.
    """
    # first import the file as table
    import_stringdb_file(ppifile, sql_conn, PAPPI_PPI_RAW_TABLE_NAME)
    
    # import ENSG<->ENSP matching
    matching.import_p2g_file(ensg2ensp_file, sql_conn, PAPPI_ENSP2ENSG_TABLE_NAME)
    
    # then create/fill the other tables (scoring, mean tissue, summary)
    init_stringdb_ppi(sql_conn)
