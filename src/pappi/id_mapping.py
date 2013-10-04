'''
Import operations for ID mapping tables. (Ensembl ENSG<->ENSP and Ensembl ENSG <-> EntreZ Gene)

@author: Patrick Flick
'''

import csv

from config import PAPPI_SQL_HGNC_FILTER_SCRIPT
from config import PAPPI_SQL_ENSEMBL2HGNC_FILTER_SCRIPT
from config import PAPPI_SQL_ENTREZ2HGNC_FILTER_SCRIPT
from config import PAPPI_SQL_UNIPROT2HGNC_FILTER_SCRIPT
import sql


COL_ENSEMBL_ID = "ensembl"
COL_ENTREZ = "entrez"
COL_HGNC_SYMB = "hgnc"
COL_HGNC_ID = "hgnc_id"
COL_UNIPROT = "uniprot"
MAPPED_IDS = [COL_ENSEMBL_ID, COL_ENTREZ, COL_HGNC_SYMB, COL_UNIPROT]
HGNC_MAPPING_TABLE_NAME = 'hgnc'
BIOMART_MAPPING_TABLE_NAME = 'biomart'


PAPPI_ENSP2ENSG_TABLE_NAME = 'ensg_to_ensp'


def import_ENSG2ENSP_file(infile, sql_conn, table=PAPPI_ENSP2ENSG_TABLE_NAME):
    """
    Imports the Ensembl protein to gene mapping ENSG<->ENSP.

    Get the ensembl ID mapping from
    http://www.ensembl.org/biomart/martview/
    using homo sapiens version 71 database and output the attributes:
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


def import_biomart_file(filename, sql_conn, table=BIOMART_MAPPING_TABLE_NAME):
    """
    Imports the mapping table from biomart.

    Go to: http://www.ensembl.org/biomart/martview

    Choose "Ensembl Genes 71" and table "Homo sapiens genes"

    Include following fields for the table:
        Ensembl Gene ID
        Associated Gene Name
        UniProt/SwissProt ID
        HGNC ID(s)
        EntrezGene ID

    Export the table as CSV (and choose "Unique results only")

    @param filename: The name of the BioMart file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the data is to be imported.
    """
    # import the CSV input as table
    cols = [COL_ENSEMBL_ID, COL_HGNC_SYMB, COL_UNIPROT, COL_HGNC_ID, COL_ENTREZ]
    sql.import_csv(filename, table + '_raw', ',', True, column_names=cols,
                   column_types=["varchar(16)"]*5, csv_quoting=None,
                   sql_conn=sql_conn)

    # remove the _HUMAN postfix from the Uniprot ID to get BioMart into the
    # same format as the HGNC data
    cur = sql_conn.cursor()
    cur.execute('DROP TABLE IF EXISTS ' + table)
    cur.execute('CREATE TABLE ' + table + ' AS '
                'SELECT '
                '' + COL_HGNC_ID + ', '
                '' + COL_HGNC_SYMB + ', '
                '' + COL_ENSEMBL_ID + ', '
                '' + COL_ENTREZ + ', '
                # remove the postfix "_HUMAN" from the Uniprot identifier
                'replace(' + COL_UNIPROT + ',"_HUMAN","") AS '
                '' + COL_UNIPROT + ' '
                'FROM ' + table + '_raw')


def import_hgnc_file(filename, sql_conn, table=HGNC_MAPPING_TABLE_NAME):
    """
    Imports the HGNC gene names table for protein-coding genes for mapping
    different gene identifiers.

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
    http://www.genenames.org/cgi-bin/hgnc_downloads?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=md_eg_id&col=md_prot_id&col=md_ensembl_id&status=Approved&status_opt=2&where=%28%28gd_pub_chrom_map+not+like+%27%25patch%25%27+and+gd_pub_chrom_map+not+like+%27%25ALT_REF%25%27%29+or+gd_pub_chrom_map+IS+NULL%29+and+gd_locus_type+%3D+%27gene+with+protein+product%27&order_by=gd_hgnc_id&format=text&limit=&submit=submit

    This should return all 19060 genes.

    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the data is to be imported.
    """

    # import the CSV file into the SQL database
    cols = [COL_HGNC_ID, COL_HGNC_SYMB, "hgnc_name", "status", COL_ENTREZ,
            COL_ENSEMBL_ID, COL_ENTREZ+"_2", COL_UNIPROT, COL_ENSEMBL_ID+"_2"]
    col_types = ["varchar(16)"]*9
    col_types[2] = "varchar(256)"
    sql.import_csv(filename, table + "_raw", '\t', True, cols, col_types,
                   None, sql_conn)

    # bring into "normal" form (i.e merge the two providers for ensembl and
    # entrez
    cur = sql_conn.cursor()

    cur.execute('DROP TABLE IF EXISTS "' + table + '"')
    query = ('CREATE TABLE ' + table + ' AS '
             'SELECT '
             '' + COL_HGNC_ID + ', '
             '' + COL_HGNC_SYMB + ', '
             # prefer the Ensembl supplied Ensembl ID
             'CASE WHEN ' + COL_ENSEMBL_ID + '_2 = "" '
             '  THEN ' + COL_ENSEMBL_ID + ' '
             '  ELSE ' + COL_ENSEMBL_ID + '_2 '
             'END AS ' + COL_ENSEMBL_ID + ', '
             # prefer the NCBI supplied Entrez ID
             'CASE WHEN ' + COL_ENTREZ + '_2 = "" '
             '  THEN ' + COL_ENTREZ + ' '
             '  ELSE ' + COL_ENTREZ + '_2 '
             'END AS ' + COL_ENTREZ + ', '
             '' + COL_UNIPROT + ' '
             'FROM ' + table + '_raw')
    cur.execute(query)
    cur.close()
    sql_conn.commit()


def create_mapping_table(from_id, to_id, sql_conn, verbose=False):
    """
    Creates the table `from_id`_2_`to_id` as an Identifier mapping table.
    Each unique identifier on the `from_id` side maps only to one identifier
    on the `to_id` side.
    """
    # check parameters, both IDs must be in the supported field IDs
    if not from_id in MAPPED_IDS and not to_id in MAPPED_IDS:
        raise ValueError("from_id and to_id must be in " + str(MAPPED_IDS))

    # name for the result table:
    table_name = from_id + "_2_" + to_id

    if verbose:
        print("Creating mapping table for " + from_id + " -> " + to_id + ":")

    # if either `to` or `from` is HGNC, then use HGNC as primary table
    if to_id == COL_HGNC_SYMB or from_id == COL_HGNC_SYMB:
        tables = ["hgnc", "biomart"]
    else:
        tables = ["biomart", "hgnc"]

    #tables = ["hgnc"]
    table1 = tables[0]

    cur = sql_conn.cursor()

    if verbose:
        print("    Creating pairwise table from: " + table1)

    # create first simple pairwise mapping table
    cur.execute('DROP TABLE IF EXISTS ' + table_name + "_pairwise")
    cur.execute('CREATE TABLE ' + table_name + '_pairwise AS '
                'SELECT ' + from_id + ', ' + to_id + ' '
                'FROM ' + table1 + ' '
                'WHERE ' + from_id + '!="" AND ' + to_id + '!=""')

    # insert pairwise matchings for those identifiers
    # in the other mapping tables, that have not yet been matched
    for table in tables[1:]:

        if verbose:
            print("    Inserting matchings from " + table)

        # create table of identifiers already matched
        cur.execute('CREATE TABLE ' + table_name + '_fromids AS '
                    'SELECT DISTINCT ' + from_id + ' '
                    'FROM ' + table_name + '_pairwise')
        # insert pairwise matching for identifers that have
        # not been matched before
        cur.execute('INSERT INTO ' + table_name + '_pairwise '
                    'SELECT ' + from_id + ', ' + to_id + ' '
                    'FROM ' + table + ' '
                    'WHERE ' + from_id + ' != "" AND ' + to_id + ' != "" '
                    'AND ' + from_id + ' NOT IN '
                    '(SELECT * FROM ' + table_name + '_fromids)')
        cur.execute('DROP TABLE ' + table_name + '_fromids')

    # remove duplicates
    # - this is done somewhat naively (choosing the identifier with MIN)
    # TODO maybe do this better?
    if verbose:
        print("    Removing duplicate and ambiguous matchings")

    cur.execute('DROP TABLE IF EXISTS ' + table_name + '_single')
    cur.execute('CREATE TABLE ' + table_name + '_single AS '
                'SELECT MIN(' + to_id + ') AS ' + to_id + ', ' + from_id + ' '
                'FROM ' + table_name + '_pairwise '
                'GROUP BY ' + from_id)

    # get statistic for #matched/#unique-ids
    # TODO get for `to_id` how many of the mismatches are because of
    #      duplicates & ambiguities removed
    if verbose:
        for i in [from_id, to_id]:
            # get the number of mapped and total gene identifiers per id
            # create table of identifiers already matched
            cur.execute('CREATE TABLE ' + table_name + '_' + i + ' AS '
                        'SELECT DISTINCT ' + i + ' '
                        'FROM ' + table_name + '_pairwise')

            # get all distinct identifier from all tables
            cur.execute('DROP TABLE IF EXISTS ' + i + '_all_ids')
            union_of_ids = " UNION ".join('SELECT ' + i + ' FROM ' + t
                                          + ' WHERE ' + i + ' != "" '
                                          for t in tables)
            cur.execute('CREATE TABLE ' + i + '_all_ids AS '
                        'SELECT DISTINCT ' + i + ' FROM '
                        '(' + union_of_ids + ')')

            # create table of unmatched identifiers
            cur.execute('CREATE TABLE ' + table_name + '_unmatched_' + i + ' '
                        'AS SELECT * FROM ' + i + '_all_ids '
                        'WHERE ' + i + ' NOT IN '
                        '( SELECT * FROM ' + table_name + '_' + i + ')')

            # get count of all
            cur.execute('SELECT COUNT() FROM ' + i + '_all_ids')
            count_from_all = cur.fetchone()[0]

            # get count of matched
            cur.execute('SELECT COUNT() FROM ' + table_name + '_' + i)
            count_from_matched = cur.fetchone()[0]

            cur.execute('DROP TABLE ' + table_name + '_' + i)

            print("    " + i + " identifiers matched:")
            print("        Matched identifiers: "
                  + str(count_from_matched) + "/" + str(count_from_all))
            print("        Unmatched identfiers: "
                  + str(count_from_all - count_from_matched))
            print("        Unmatched identifiers are now saved in the `"
                  + table_name + "_unmatched_" + i + "` SQL table.")

    # close cursor and commit changes to SQL server
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



def import_mappings(hgnc_file, biomart_file, sql_conn):
    """
    @param hgnc_file:   The HGNC mapping file, see the documentation of
                        the function import_hgnc_file() in this module.
    @param biomart_file:    The Biomart mapping file, see the documentation
                            of the import_biomart_file() function in this
                            module.
    @param sql_conn: The SQL connection to be used.
    """
    
    # first import the hgnc file
    import_hgnc_file(hgnc_file, sql_conn)
    
    # then import the biomart file
    import_biomart_file(biomart_file, sql_conn)
    
    # then initialize all the mappings (creating mapping tables to hgnc from other IDs)
    
    # Ensembl -> HGNC 
    execute_script(PAPPI_SQL_ENSEMBL2HGNC_FILTER_SCRIPT, sql_conn)
    
    # Uniprot -> HGNC
    execute_script(PAPPI_SQL_UNIPROT2HGNC_FILTER_SCRIPT, sql_conn)
    
    # Entrez -> HGNC 
    execute_script(PAPPI_SQL_ENTREZ2HGNC_FILTER_SCRIPT, sql_conn)


def import_hgnc_entrez2ensembl(infile, sql_conn):
    """
    Imports the HGNC table for protein-coding genes for mapping of Entrez Gene ID to
    Ensembl Gene ID, then filters and merges the data to create an Entrez to Ensembl mapping
    table.
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    """
    # first import the HGNC file
    import_hgnc_entrez2ensembl_file(infile, sql_conn)
    
    # then filter and create the ID mapping table
    init_entrez2ensembl(sql_conn)
