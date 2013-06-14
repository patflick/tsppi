import sqlite3
import csv


filename= "../../hpa/data/tissues/intersected/cerebral.cortex_neuronal.cells.csv"


con = sqlite3.Connection('newdb.sqlite')
cur = con.cursor()
# "Gene","Tissue","Cell type","Level","Expression type","Reliability"
cur.execute('DROP TABLE "cerebral";');
cur.execute('CREATE TABLE "cerebral" ("Gene" varchar(16), "Tissue" varchar(256), "Cell.type" varchar(256), "Level" varchar(32), "Expression.type" varchar(12), "Reliability" varchar(32));')

f = open(filename)
csv_reader = csv.reader(f, delimiter=',',quoting=csv.QUOTE_ALL)

cur.executemany('INSERT INTO cerebral VALUES (?, ?, ?, ?, ?, ?)', csv_reader)

cur.execute()

cur.close()
con.commit()
con.close()
f.close()