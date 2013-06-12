'''
Import and filtering operations for PPI data.

@author: Patrick Flick
'''

import csv
import re

from . import mapping
from config import PAPPI_SQL_STRINGDB_FILTER_SCRIPT
from config import PAPPI_SQL_CCSB_FILTER_SCRIPT
from config import PAPPI_SQL_MMC_FILTER_SCRIPT
from config import PAPPI_SQL_EDGE_EXPR_SCRIPT
from sql import execute_script

# TODO put all constants in an own module
PAPPI_STRINGDB_RAW_TABLE_NAME = 'stringdb_raw'
PAPPI_CCSB_RAW_TABLE_NAME = 'ccsb_raw'
PAPPI_MMC_RAW_TABLE_NAME = 'mmc_raw'

PAPPI_ENSP2ENSG_TABLE_NAME = 'ensg_to_ensp'


#################################
# String-DB import and convert
#################################


def import_stringdb_file(infile, sql_conn, table=PAPPI_STRINGDB_RAW_TABLE_NAME):
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
    execute_script(PAPPI_SQL_STRINGDB_FILTER_SCRIPT, sql_conn)



def import_stringdb(ppifile, ensg2ensp_file, sql_conn):
    """
    Imports and initializes the string-db PPI network. Also the given ENSG<->ENSP mapping file
    is imported, which is used during the string-db initialization.
    
    @param ppifile: The file handle for the PPI network file, downloaded from string-db.org
    @param ensg2ensp_file: The file handle for the ENSG<->ENSP mapping to be imported.
    @param sql_conn: The SQL connection to be used.
    """
    # first import the file as table
    import_stringdb_file(ppifile, sql_conn, PAPPI_STRINGDB_RAW_TABLE_NAME)
    
    # import ENSG<->ENSP mapping
    mapping.import_ENSG2ENSP_file(ensg2ensp_file, sql_conn, PAPPI_ENSP2ENSG_TABLE_NAME)
    
    # then perform the ENSG<->ENSP mapping, returning the PPI network
    # that uses ENSG IDs (this is the same format that results from the CCSB import & init functions)
    init_stringdb_ppi(sql_conn)


#################################
# CCSB import and convert
#################################

def import_ccsb_file(infile, sql_conn, table=PAPPI_CCSB_RAW_TABLE_NAME):
    """
    Imports the CCSB PPI network from the tab separated file. The PPI is imported
    using the `sql_conn` SQL connection into a new table given by `table`.
    
    This imports the data from:
        http://interactome.dfci.harvard.edu/H_sapiens/index.php?page=newrelease
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the CCSB data is to be imported.
    """
    # initialize the cursor object
    cur = sql_conn.cursor()
    
    # create table for the raw HPA data:
    # GENE_IDA    SYMBOL_A    GENE_IDB    SYMBOL_B
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("Gene_IDA" int, "Symbol_A" varchar(16), "Gene_IDB" int, "Symbol_B" varchar(16))')
    
    # get csv reader for the CCSB file
    csv_reader = csv.reader(infile, delimiter='\t',quoting=csv.QUOTE_NONE)
    
    # ignore header line
    csv_reader.next()
    
    # insert all lines
    cur.executemany('INSERT INTO "' + table + '" VALUES (?, ?, ?, ?)', csv_reader)

    # close cursor and commit
    cur.close()
    sql_conn.commit()
    

def init_ccsb_ppi(sql_conn):
    """
    Initializes the imported CCSB PPI data. This means the Entrez IDs
    are mapped to Ensembl IDs using the `entrez_to_ensembl` table.
    
    @param sql_conn: The SQL connection to be used.
    """
    execute_script(PAPPI_SQL_CCSB_FILTER_SCRIPT, sql_conn)



def import_ccsb(ccsb_file, hgnc_file, sql_conn):
    """
    Imports and initializes the CCSB PPI network. Also the given HGNC mapping file
    is imported, which is used for mapping the Entrez Gene IDs of the CCSB data
    to Ensembl Gene IDs used in HPA.
    
    @param ppifile: The file handle for the PPI network file, downloaded from CCSB
    @param ensg2ensp_file: The file handle for the ENSG<->ENSP mapping to be imported.
    @param sql_conn: The SQL connection to be used.
    """
    # first import the file as table
    import_ccsb_file(ccsb_file, sql_conn, PAPPI_CCSB_RAW_TABLE_NAME)
    
    # then perform the entrez to ensembl mapping, returning the PPI network
    # that uses ENSG IDs (this is the same format that results from the string-db import & init functions)
    init_ccsb_ppi(sql_conn)


#################################
# Cell: human complex PPI
#################################

def import_mmc_file(infile, sql_conn, table=PAPPI_MMC_RAW_TABLE_NAME):
    """
    Imports the MMC (A Census of Human Soluble Protein Complexes, Havugimana et al.)
     PPI network from the tab separated file. The PPI is imported
    using the `sql_conn` SQL connection into a new table given by `table`.
    
    This imports the data from:
        http://interactome.dfci.harvard.edu/H_sapiens/index.php?page=newrelease
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the CCSB data is to be imported.
    """
    # initialize the cursor object
    cur = sql_conn.cursor()
    
    # create table for the raw HPA data:
    # GENE_IDA    SYMBOL_A    GENE_IDB    SYMBOL_B
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("Gene1" varchar(10), "Gene2" varchar(10))')
    
    # get csv reader for the CCSB file
    csv_reader = csv.reader(infile, delimiter='\t',quoting=csv.QUOTE_NONE)
    
    # ignore header line
    csv_reader.next()
    
    # insert all lines
    cur.executemany('INSERT INTO "' + table + '" VALUES (?, ?)', csv_reader)

    # close cursor and commit
    cur.close()
    sql_conn.commit()



def import_mmc(mmc_file, sql_conn):
    """
    Imports and initializes the CCSB PPI network. Also the given HGNC mapping file
    is imported, which is used for mapping the Entrez Gene IDs of the CCSB data
    to Ensembl Gene IDs used in HPA.
    
    @param ppifile: The file handle for the PPI network file, downloaded from CCSB
    @param ensg2ensp_file: The file handle for the ENSG<->ENSP mapping to be imported.
    @param sql_conn: The SQL connection to be used.
    """
    # first import the file as table
    import_mmc_file(mmc_file, sql_conn, PAPPI_MMC_RAW_TABLE_NAME)
    
    # run script that creates the PPI with HGNC symbols, rather than UniprotIDs
    execute_script(PAPPI_SQL_MMC_FILTER_SCRIPT, sql_conn)



#################################
# Processing common to all PPIs
#################################

def init_edge_expression(sql_conn):
    """
    Creates tissue_expr and edge_expression tables for the current ppi.
    
    @param sql_conn: The SQL connection to be used.
    """
    execute_script(PAPPI_SQL_EDGE_EXPR_SCRIPT, sql_conn)
