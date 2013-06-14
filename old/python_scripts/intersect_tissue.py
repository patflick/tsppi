

import csv
import sqlite3
import argparse

DATABASE_NAME='hpa_tisses.sqlite'
MERGE_SQL_SCRIPT='merge_tissues.sql'
LEVEL_MAPPING_TABLE='LevelMapping'


def sql_create_mapping(database, table=LEVEL_MAPPING_TABLE):
    con = sqlite3.Connection(database)
    cur = con.cursor()
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("Level" varchar(16), "Score" int);')
    cur.execute('INSERT INTO "' + table + '" VALUES ("High", 3);')
    cur.execute('INSERT INTO "' + table + '" VALUES ("Medium", 2);')
    cur.execute('INSERT INTO "' + table + '" VALUES ("Low", 1);')
    cur.execute('INSERT INTO "' + table + '" VALUES ("None", 0);')
    cur.close()
    con.commit()
    con.close()

def sql_import_hpa_tissue(infile, database, table):
    con = sqlite3.Connection(database)
    cur = con.cursor()
    # "Gene","Tissue","Cell type","Level","Expression type","Reliability"
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("Gene" varchar(16), "Tissue" varchar(256), "Cell.type" varchar(256), "Level" varchar(32), "Expression.type" varchar(12), "Reliability" varchar(32));')
    
    
    csv_reader = csv.reader(infile, delimiter=',',quoting=csv.QUOTE_ALL)
    header = csv_reader.next()
    
    cur.executemany('INSERT INTO "' + table + '" VALUES (?, ?, ?, ?, ?, ?)', csv_reader)
    
    cur.close()
    con.commit()
    con.close()

def sql_merge(database):
    f = open(MERGE_SQL_SCRIPT, 'r')
    qry = f.read()
    f.close()
    # open SQL connection
    con = sqlite3.Connection(database)
    cur = con.cursor()
    # execute query
    cur.execute('DROP TABLE IF EXISTS "merged_tissues"')
    cur.execute(qry)
    # commit and close SQL
    cur.close()
    con.commit()
    con.close()

def sql_dump_csv(database, table, outfile):
    # open SQL connection
    con = sqlite3.Connection(database)
    cur = con.cursor()
    
    cols = []
    for row in cur.execute('PRAGMA table_info("' + table + '")'):
        cols.append(row[1])
        
    print cols
    
    wr = csv.writer(outfile, delimiter=' ', quoting=csv.QUOTE_NONE)
    wr.writerow(cols)
    
    for row in cur.execute('SELECT * FROM "' + table + '"'):
        wr.writerow(row)
    
    cur.close()
    con.commit()
    con.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Intersects two HPA tissue data sets,\nand outputs only those genes that are scored with `APE`\nin both datasets.")
    parser.add_argument('inputfile1',type=argparse.FileType("r"))
    parser.add_argument('inputfile2',type=argparse.FileType("r"))
    parser.add_argument('outputfile1',type=argparse.FileType("w"))
    #parser.add_argument('outputfile2',type=argparse.FileType("w"))
    #parser.add_argument('-s','--min-score', dest='min_score',help='the minimum APE reliability score. genes that score less than this in either one of the sets, are not included in the intersection. (default=Medium)', default='Medium')
    #parser.add_argument('-o','--output-genes', dest='out_genes',default=None,type=argparse.FileType('w'),help='outputs the genes that are in the intersection, separated by newline')
    args = parser.parse_args()
    
    sql_create_mapping(DATABASE_NAME)
    sql_import_hpa_tissue(args.inputfile1, DATABASE_NAME, 'tissue1')
    sql_import_hpa_tissue(args.inputfile2, DATABASE_NAME, 'tissue2')
    sql_merge(DATABASE_NAME)
    sql_dump_csv(DATABASE_NAME, 'merged_tissues', args.outputfile1)
    