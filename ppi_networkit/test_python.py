
import sys

sys.path.append("./cython")
#import NetworKit
import ppi_networkit

sqlio = ppi_networkit.SQLiteIO("/home/patrick/dev/bio/data/test_matching.sqlite")

tg1 = sqlio.load_tsppi_graph("ccsb", "hpa")

name = tg1.getPpiName()
