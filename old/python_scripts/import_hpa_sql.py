import sqlite3
import csv

filename= "../../hpa/data/normal_tissue.csv"
table_name="hpa_normal_tissue"

con = sqlite3.Connection('hpaDB.sqlite')
cur = con.cursor()
# "Gene","Tissue","Cell type","Level","Expression type","Reliability"
cur.execute('DROP TABLE IF EXISTS "' + table_name + '";');
cur.execute('CREATE TABLE "' + table_name + '" ("Gene" varchar(16), "Tissue" varchar(256), "Cell.type" varchar(256), "Level" varchar(32), "Expression.type" varchar(12), "Reliability" varchar(32));')

f = open(filename)
csv_reader = csv.reader(f, delimiter=',',quoting=csv.QUOTE_ALL)

cur.executemany('INSERT INTO "' + table_name + '" VALUES (?, ?, ?, ?, ?, ?)', csv_reader)

cur.close()
con.commit()
con.close()
f.close()