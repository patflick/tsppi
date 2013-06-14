


import csv
import sqlite3
import argparse
import math
import networkx as nx

DATABASE_NAME='hpa_tisses.sqlite'
LEVEL_MAPPING_TABLE='LevelMapping'
PPI_TABLE='ppi'
# TODO put this as argument
PPI_FILE='../../hpa/data/tissues/intersected/ppi.csv'
# TODO as argument, clean up code!
DIFF_OUTPUT_CSV='../../hpa/data/tissues/intersected/diff_score.csv'

# TODO put this function somewhere else
def sql_import_ppi(infile, database, table):
    con = sqlite3.Connection(database)
    cur = con.cursor()
    # "Gene","Tissue","Cell type","Level","Expression type","Reliability"
    cur.execute('DROP TABLE IF EXISTS "' + table + '";');
    cur.execute('CREATE TABLE "' + table + '" ("Gene1" varchar(16), "Gene2" varchar(16), "Score" int);')
    
    
    csv_reader = csv.reader(infile, delimiter=' ',quoting=csv.QUOTE_NONE)
    header = csv_reader.next()
    
    cur.executemany('INSERT INTO "' + table + '" VALUES (?, ?, ?)', csv_reader)
    cur.close()
    con.commit()
    con.close()

def create_graph(database, table):
    con = sqlite3.Connection(database)
    cur = con.cursor()
    G = nx.Graph()
    # "Gene","Tissue","Cell type","Level","Expression type","Reliability"
    print "Getting gene nodes ..."
    query = '''
        SELECT
            a.Gene AS Gene,
            b.Score as Score
        FROM "%s" AS a
        INNER JOIN "%s" as b
            ON a.Level = b.Level
        WHERE
            (
                a.Reliability = "High"
            OR
                a.Reliability = "Medium"
            )
            AND
                a."Expression.type" = "APE"
        ''' % (table, LEVEL_MAPPING_TABLE)
    
    for row in cur.execute(query):
        G.add_node(row[0])
        G[row[0]]['score'] = row[1]
        
    print "Getting ppi edges ..."
    query = '''
        SELECT
            Gene1,
            Gene2
        FROM ppi
    '''
    
    for row in cur.execute(query):
        G.add_edge(row[0], row[1])
        
    print "graph has # nodes= " + str(G.number_of_nodes())
    print "graph has # edges= " + str(G.number_of_edges())
    print ""
    print "Calculating neighbor-score..."
    
    for node, edges in G.adjacency_iter():
        node_sum = 0;
        node_num = 0;
        for node_neighbor, edge_atr in edges.items():
            if node_neighbor != 'score':
                node_sum += G[node_neighbor]['score']
                node_num += 1
        node_sum += G[node]['score']
        node_num += 1
        G[node]['newscore'] =  1.0* node_sum / math.sqrt(1.0*node_num)
    
    
    newtable = table + "_scored"
    print "Inserting into new table %s ..." % newtable
    cur.execute('DROP TABLE IF EXISTS "%s"' % newtable)
    cur.execute('CREATE TABLE "%s" ("Gene" varchar(16), "Score" float);' % newtable)
    
    for node, edges in G.adjacency_iter():
        cur.execute('INSERT INTO "%s" VALUES ("%s", %f);' % (newtable, node, G[node]['newscore']))
        #print 'INSERT INTO "%s" VALUES (Gene="%s",Score=%f)' % (newtable, node, G[node]['newscore'])
    
    cur.close()
    con.commit()
    con.close()
    
# TODO this file is double, also in intersect_tissue.py 
# TODO clean up, module-arize everything!
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
    
def sql_create_diff_score_table(database, table1, table2, newtable):
    con = sqlite3.Connection(database)
    cur = con.cursor()
    
    cur.execute('DROP TABLE IF EXISTS "%s"' % newtable)
    
    query = '''
        CREATE TABLE "%s" AS
        SELECT
            a.Gene,
            a.Score AS Score1,
            b.Score AS Score2,
            a.Score - b.Score AS DiffScore
        FROM
            "%s" AS a
        INNER JOIN
            "%s" AS b
        ON a.Gene = b.Gene
    ''' % (newtable, table1, table2)
    
    cur.execute(query)
    
    cur.close()
    con.commit()
    con.close()
    
if __name__ == "__main__":
    
    print "importing ppi into DB..."
    ppifile = open(PPI_FILE, 'r')
    sql_import_ppi(ppifile, DATABASE_NAME, PPI_TABLE)
    ppifile.close()
    
    print "========================"
    print " Scoring tissue 1"
    print "========================"
    create_graph(DATABASE_NAME, 'tissue1')
    print "========================"
    print " Scoring tissue 2"
    print "========================"
    create_graph(DATABASE_NAME, 'tissue2')
    print "========================"
    print " Creating diff table"
    print "========================"
    sql_create_diff_score_table(DATABASE_NAME, 'tissue1_scored', 'tissue2_scored', 'DiffScore')
    print "========================"
    print " Dumping diff table"
    print "========================"
    outfile = open(DIFF_OUTPUT_CSV, 'w')
    sql_dump_csv(DATABASE_NAME, 'DiffScore', outfile)
    outfile.close()
    
    
    
    