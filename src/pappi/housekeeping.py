

import csv

def import_entrez_file(infile, sql_conn, table):
    """
    imports the housekeeping genes as entrez (rather than ensembl) ids,
    which need conversion afterwards (init function)
    
    @param infile: The opened file handle of the file to be imported.
    @param sql_conn: The SQL connection to be used.
    @param table: The SQL table name into which the CCSB data is to be imported.
    """
    # initialize the cursor object
    cur = sql_conn.cursor()
    
    # create table for the raw HPA data:
    # GENE_IDA    SYMBOL_A    GENE_IDB    SYMBOL_B
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("Gene_ID" int)')
    
    # get csv reader for the CCSB file
    csv_reader = csv.reader(infile, delimiter='\t',quoting=csv.QUOTE_NONE)
    
    #print csv_reader
    #for item in csv_reader:
    #    print item

    # insert all lines
    cur.executemany('INSERT INTO "' + table + '" VALUES (?)', csv_reader)

    # close cursor and commit
    cur.close()
    sql_conn.commit()
    

def translate_entrez_2_ensembl(sql_conn, src_table, dst_table):
    """
    takes the src_table which is a list of Entrez Gene ids and
    creates dst_table with the corresponding Ensembl ids

    @param sql_conn: The SQL connection to be used.
    """

    cur = sql_conn.cursor()

    cur.execute("DELETE FROM entrez_to_ensembl WHERE EntrezID = '' OR EnsemblID = ''")
    cur.execute("DROP TABLE IF EXISTS " + dst_table)
    cur.execute("CREATE TABLE " + dst_table + " AS SELECT a.EnsemblID as Gene FROM entrez_to_ensembl as a INNER JOIN " + src_table + " AS b ON a.EntrezID = b.Gene_ID")

    cur.close()
    sql_conn.commit()

