'''
Import operations for ID mapping tables. (Ensembl ENSG<->ENSP and Ensembl ENSG
<-> EntreZ Gene)

@author: Patrick Flick
'''
# TODO list for this file:
#  - ensg <-> ensp mapping directly via the biomart table
#  - (maybe) replace all the DROP TABLE IF EXISTS; CREATE TABLE
#            with CREATE TABLE IF NOT EXISTS
#  - find rationale for MIN(id) to resolve ambiguties
#  - idk what else yet

import csv

from . import sql


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
    cur.execute('DROP TABLE IF EXISTS "' + table + '";')
    cur.execute('CREATE TABLE "' + table + '" ("ENSG" varchar(16), "ENSP" varchar(256));')

    # get csv reader for the hpa file
    csv_reader = csv.reader(infile, delimiter=',', quoting=csv.QUOTE_ALL)
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


def map_identifier(from_table, from_cols, from_id, to_table, to_id,
                   sql_conn, verbose=False):
    """
    Takes any SQL table and creates a new SQL table, replacing the given column
    of gene identifiers by a new column of gene identifiers of another kind.
    E.g. Mapping all Ensembl identifiers to HGNC Symbols.

    @param from_table:  The original SQL table.
    @param from_cols:   The columns to be mapped to the new identifier.
    @param from_id:     The identifier type of the original table.
                        (E.g. 'hgnc', 'ensembl', 'uniprot', 'entrez')
    @param to_table:    The name for the new SQL table to be created.
    @param to_id:       The identfier type to be mapped to.
                        (E.g. 'hgnc', 'ensembl', 'uniprot', 'entrez')
    @param sql_conn:    The SQL connection object.
    @param verbose:     Whether or not to print out debug information, this
                        may result in an elongated run time (due to
                        calculating statistics). Default: False
    """
    if verbose:
        print("Matching Table '" + from_table + "' -> '" + to_table + "' "
              "by mapping identifiers " + from_id + " -> " + to_id)

    # check if the needed mapping table is already present
    # otherwise first create it
    mapping_table = from_id + '_2_' + to_id
    if not sql.table_exists(mapping_table, sql_conn):
        if verbose:
            print("    Mapping table for " + from_id + " -> " + to_id
                  + " doesn't yet exists.")
        # create new mapping table
        create_mapping_table(from_id, to_id, sql_conn, verbose)

    # get all the column names in order to be able to generate the
    # SQL statement
    col_names = sql.get_column_names(from_table, sql_conn)

    # get a SQl cursor object
    cur = sql_conn.cursor()

    # get data type of the from_id of the mapping table
    # for data conversions (if the according column in the from_table is
    # not converted to the same data type, then the indexes are not used,
    # resulting in very very poor performance)
    cur.execute('SELECT TYPEOF(' + from_id + ') FROM ' + mapping_table)
    from_id_type = cur.fetchone()[0]

    # construct the fields for the result table and
    # the JOIN statement composed of INNER JOINs
    mapped_col_count = 0
    result_fields = []
    join_commands = []
    for col in col_names:
        if col in from_cols:
            # create map alias in the form of map1, map2, map3, ...
            # for the join statements and the correct selection
            # of the field
            mapped_col_count += 1
            map_alias = 'map' + str(mapped_col_count)
            # construct the result field statement (e.g. map1.hgnc AS gene1)
            result_fields.append(map_alias + '.' + to_id + ' AS ' + col)
            # construct the inner join statement
            # (e.g. INNER JOIN entrez_2_hgnc AS map1 ON map1.entrez = src.gene1)
            join_commands.append('INNER JOIN ' + mapping_table + ' AS '
                                 + map_alias + ' ON ' + map_alias + '.'
                                 # CAST the col to the type of the mapping
                                 # table (otherwise INDEX is not used)
                                 + from_id + ' = ' + 'CAST(src.' + col + ' AS '
                                 + from_id_type + ')')
        else:
            # construct result field stmt (e.g. src.descr AS descr)
            result_fields.append('src.' + col + ' AS ' + col)

    # contruct the SQL query to do the inner joins
    sqlquery = ('SELECT '
                + ", ".join(result_fields) + ' '
                'FROM ' + from_table + ' AS src '
                + " ".join(join_commands))

    # for performance reasons: create an index for the from_id column of the
    # mapping table (in case it does not yet exist)
    cur.execute('CREATE INDEX IF NOT EXISTS ' + mapping_table + '_' + from_id
                + '_index ON ' + mapping_table + '(' + from_id + ')')

    # run the inner join
    cur.execute('DROP TABLE IF EXISTS ' + to_table)
    cur.execute('CREATE TABLE ' + to_table + ' AS ' + sqlquery)

    # verbose debug output: get table and statistics of non-matched rows
    if verbose:
        # get all from_ids from all used columns
        union_of_ids = " UNION ".join('SELECT ' + c + ' AS id FROM '
                                      + from_table
                                      for c in from_cols)
        cur.execute('DROP TABLE IF EXISTS ' + from_table + '_orig_ids')
        cur.execute('CREATE TABLE ' + from_table + '_orig_ids AS '
                    'SELECT DISTINCT id FROM (' + union_of_ids + ')')

        # create table of non-matched ids
        cur.execute('DROP TABLE IF EXISTS ' + to_table + '_unmatched_ids')
        cur.execute('CREATE TABLE ' + to_table + '_unmatched_ids AS '
                    'SELECT DISTINCT id FROM ' + from_table + '_orig_ids '
                    'WHERE id NOT IN ('
                    '   SELECT ' + from_id + ' FROM ' + mapping_table + ')')

        # create table made of rows which did not map
        where_or = ") OR (".join(c + ' IN (SELECT id FROM ' + to_table
                                 + '_unmatched_ids)' for c in from_cols)
        where_cond = "(" + where_or + ")"
        cur.execute('DROP TABLE IF EXISTS ' + to_table + '_unmatched_rows')
        cur.execute('CREATE TABLE ' + to_table + '_unmatched_rows AS '
                    'SELECT * FROM ' + from_table + ' '
                    'WHERE ' + where_cond)

        # get statistics (number of matched and unmatched rows and ids)
        cur.execute('SELECT COUNT(*) FROM ' + from_table)
        nrows_from = cur.fetchone()[0]

        cur.execute('SELECT COUNT(*) FROM ' + to_table)
        nrows_to = cur.fetchone()[0]

        cur.execute('SELECT COUNT(*) FROM ' + to_table + '_unmatched_rows')
        nrows_unmatched = cur.fetchone()[0]

        cur.execute('SELECT COUNT(*) FROM ' + from_table + '_orig_ids')
        num_from_ids = cur.fetchone()[0]

        cur.execute('SELECT COUNT(*) FROM ' + to_table + '_unmatched_ids')
        num_unmatched_ids = cur.fetchone()[0]

        if (nrows_from - nrows_to != nrows_unmatched):
            raise RuntimeError("Implementation Error - numbers don't match up.")

        # print out all the statistics
        print("Matched Table '" + from_table + "' -> '" + to_table + "'")
        print("    Matched " + str(num_from_ids) + " " + from_id + " to "
              + str(num_from_ids - num_unmatched_ids) + " " + to_id
              + " identifiers (" + str(num_unmatched_ids) + " unmatched)")
        print("    Successfully matched rows: " + str(nrows_to) + "/"
              + str(nrows_from))
        if (num_unmatched_ids > 0):
            print("    Unmatched rows ( = " + str(nrows_unmatched)
                  + " ) are now available in the SQL table `" + to_table
                  + "_unmatched_rows`")
            print("    Unmatched identfiers ( = " + str(num_unmatched_ids)
                  + " ) are now available in the SQL table `" + to_table
                  + "_unmatched_ids`")

        # clean up somewhat
        cur.execute('DROP TABLE ' + from_table + '_orig_ids')

    # close cursor and commit to server
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

    # TODO maybe check if the tables even exists and quit with an error if not

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
    # TODO how else can ambiguities be resolved?
    # TODO how does this change/influence the PPIs for example?
    if verbose:
        print("    Removing duplicate and ambiguous matchings")

    cur.execute('DROP TABLE IF EXISTS ' + table_name + '_single')
    cur.execute('CREATE TABLE ' + table_name + '_single AS '
                'SELECT MIN(' + to_id + ') AS ' + to_id + ', ' + from_id + ' '
                'FROM ' + table_name + '_pairwise '
                'GROUP BY ' + from_id)

    # get statistic for #matched/#unique-ids
    if verbose:
        for i in [from_id, to_id]:
            # get the number of mapped and total gene identifiers per id
            # create table of identifiers already matched
            cur.execute('CREATE TABLE ' + table_name + '_' + i + ' AS '
                        'SELECT DISTINCT ' + i + ' '
                        'FROM ' + table_name + '_single')

            # get all distinct identifier from all tables
            cur.execute('DROP TABLE IF EXISTS ' + i + '_all_ids')
            union_of_ids = " UNION ".join('SELECT ' + i + ' FROM ' + t
                                          + ' WHERE ' + i + ' != "" '
                                          for t in tables)
            cur.execute('CREATE TABLE ' + i + '_all_ids AS '
                        'SELECT DISTINCT ' + i + ' FROM '
                        '(' + union_of_ids + ')')

            # create table of unmatched identifiers
            cur.execute('DROP TABLE IF EXISTS ' + table_name + '_unmatched_'
                        + i + ' ')
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

            # get count of non-matched because of removal of replicates
            unmatched_by_rm_dupl = -1
            if i == to_id:
                cur.execute('SELECT COUNT() FROM ('
                            'SELECT DISTINCT ' + i + ' '
                            'FROM ' + table_name + '_pairwise '
                            'WHERE ' + i + ' NOT IN ('
                            'SELECT DISTINCT ' + i + ' '
                            'FROM ' + table_name + '_single))')
                unmatched_by_rm_dupl = cur.fetchone()[0]

            cur.execute('DROP TABLE ' + table_name + '_' + i)

            print("    " + i + " identifiers matched:")
            print("        Matched identifiers: "
                  + str(count_from_matched) + "/" + str(count_from_all))
            print("        Unmatched identfiers: "
                  + str(count_from_all - count_from_matched))
            if unmatched_by_rm_dupl >= 0:
                print("            because of removal of duplicates: "
                      + str(unmatched_by_rm_dupl) + "/"
                      + str(count_from_all-count_from_matched))
            print("        Unmatched identifiers are now saved in the `"
                  + table_name + "_unmatched_" + i + "` SQL table.")

    # create the real mapping table
    cur.execute('DROP TABLE IF EXISTS ' + table_name)
    cur.execute('CREATE TABLE ' + table_name + ' AS '
                'SELECT ' + from_id + ', ' + to_id + ' '
                'FROM ' + table_name + '_single '
                'ORDER BY ' + from_id)

    # close cursor and commit changes to SQL server
    cur.close()
    sql_conn.commit()
