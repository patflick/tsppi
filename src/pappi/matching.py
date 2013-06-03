'''
Import operations for ID matching tables. (Ensembl ENSG<->ENSP and Ensembl ENSG <-> EntreZ Gene)

@author: Patrick Flick
'''

import csv

from config import PAPPI_SQL_HGNC_FILTER_SCRIPT
from sql import execute_script

PAPPI_ENSP2ENSG_TABLE_NAME = 'ensg_to_ensp'
PAPPI_HGNC_MAPPING_TABLE_NAME = 'hgnc_raw'

def import_ENSG2ENSP_file(infile, sql_conn, table=PAPPI_ENSP2ENSG_TABLE_NAME):
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
    


def import_hgnc_entrez2ensembl_file(infile, sql_conn, table=PAPPI_HGNC_MAPPING_TABLE_NAME):
    """
    Imports the HGNC table for protein-coding genes for mapping of Entrez Gene ID to
    Ensembl Gene ID.
    
    Get data from: http://www.genenames.org/cgi-bin/hgnc_stats
    Goto `Locus Group`: "protein-coding gene" and click "Custom".
    Choose only the Columns:
        - HGNC ID
        - Approved Symbol
        - Approved Name
        - Status
        - Entrez Gene ID
        - Ensembl Gene ID
        (and from external sources)
        - Entrez Gene ID (supplied by NCBI)
        - UniProt ID (supplied by UniProt)
        - Ensembl ID (supplied by Ensembl)
    Make sure to deselect (exclude) the status: "Entry and Symbol Withdrawn"
    
    Full URL to results:
    http://www.genenames.org/cgi-bin/hgnc_downloads?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=md_eg_id&col=md_prot_id&col=md_ensembl_id&status=Approved&status_opt=2&where=%28%28gd_pub_chrom_map+not+like+%27%25patch%25%27+and+gd_pub_chrom_map+not+like+%27%25ALT_REF%25%27%29+or+gd_pub_chrom_map+IS+NULL%29+and+gd_locus_type+%3D+%27gene+with+protein+product%27&order_by=gd_hgnc_id&format=text&limit=&hgnc_dbtag=on&submit=submit
        
    This should return all 19060 genes.
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the data is to be imported.
    """
    # initialize the cursor object
    cur = sql_conn.cursor()
    
    # create table for the raw HPA data:
    # HGNC ID    Approved Symbol    Approved Name    Entrez Gene ID    Ensembl Gene ID    Entrez Gene ID    Ensembl ID
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("HGNC.ID" varchar(16), "HGNC.Symbol" varchar(16), "HGNC.Name" varchar(256), "HGNC.Status", "HGNC.EntrezID" int, "HGNC.EnsemblID" varchar(16), "NCBI.EntrezID" int, "UniprotID" varchar(10), "Ensembl.EnsemblID" varchar(16));')
    
    # get csv reader for the hpa file
    csv_reader = csv.reader(infile, delimiter='\t',quoting=csv.QUOTE_NONE)
    # ignore header line
    csv_reader.next()
    
    # insert all lines
    cur.executemany('INSERT INTO "' + table + '" VALUES (?,?,?,?,?,?,?,?,?)', csv_reader)
    
    # close cursor and commit
    cur.close()
    sql_conn.commit()


def init_entrez2ensembl(sql_conn):
    """
    Filters the imported HGNC dataset and outputs a new table `entrez_to_ensembl`
    with a mapping from Entrez Gene IDs to Ensembl Gene IDs.
    
    @todo TODO There are some unmatched genes on both sides (19060 - 18895)
    
    @param sql_conn: The SQL connection to be used.
    """
    execute_script(PAPPI_SQL_HGNC_FILTER_SCRIPT, sql_conn)



def import_hgnc_entrez2ensembl(infile, sql_conn):
    """
    Imports the HGNC table for protein-coding genes for mapping of Entrez Gene ID to
    Ensembl Gene ID, then filters and merges the data to create an Entrez to Ensembl matching
    table.
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    """
    # first import the HGNC file
    import_hgnc_entrez2ensembl_file(infile, sql_conn)
    
    # then filter and create the ID matching table
    init_entrez2ensembl(sql_conn)
