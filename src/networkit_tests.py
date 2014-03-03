# this file is for interactive tests
# thus it just loads the configuration and
# database connection and enables tab completion

import os
import re
import pappi.id_mapping
import pappi.sql
from pappi.data_config import *

##############################################################
#  TODO: temporary content (needs to be put into a module)   #
##############################################################

# histogram function using counter
from collections import Counter


def cluster_hist(clusters):
    c = Counter()
    cluster_sizes = clusters.subsetSizes()
    for size in cluster_sizes:
        c[size] += 1
    return c


def print_sorted_hist(hist):
    print("size\t#clusters")
    for size, num in iter(sorted(hist.items())):
        print(str(size) + ":\t" + str(num))

import sys
sys.path.append("/home/patrick/dev/bio/goatools")
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag

def load_go_associations(sql_conn, only_genes=None):
    cur = sql_conn.cursor()
    cur.execute("SELECT * FROM go_gene_assoc")
    assoc = {}
    for row in cur.fetchall():
        gene = row[0]
        go_term = row[1]
        if only_genes and gene not in only_genes:
            continue
        if gene in assoc:
            assoc[gene].add(go_term)
        else:
            assoc[gene] = set([go_term])
    return assoc


def enrichment_scoring(cluster, all_genes, assoc):
    print("got cluster of size " + str(len(cluster)) + " and #genes: " + str(len(all_genes)))
    print("Performing gene enrichment study:")
    obo_dag = GODag(obo_file=GO_OBO_FILE)
    population = set(all_genes)
    study = set(cluster)
    ge = GOEnrichmentStudy(population, assoc, obo_dag, alpha=0.05,
            study=study, methods=["bonferroni", "sidak", "holm"])
    ge.print_summary(pval=0.05)


def score_clusters(clusters, tsppi, all_genes, assoc):
    cluster_sizes = clusters.subsetSizeMap()
    # loop through all the clusters, get the associated gene names
    # and call the GO enrichment scoring of the cluster
    for cluster, size in cluster_sizes.items():
        cluster_nodes = clusters.getMembers(cluster)
        cl_node_names = [tsppi.getGeneName(x) for x in cluster_nodes]
        enrichment_scoring(cl_node_names, all_genes, assoc)

##############################
# get database connection
##############################

# get new database connection
con = pappi.sql.get_conn(DATABASE)



#################################
#  Import PPI NetworKit module  #
#################################

#sys.path.append("/home/patrick/dev/bio/NetworKit-Flick/cython")
sys.path.append("/home/patrick/dev/bio/ppi_networkit/cython")
#import NetworKit
import ppi_networkit

sqlio = ppi_networkit.SQLiteIO(DATABASE)

tsppi = sqlio.load_tsppi_graph("string", "gene_atlas")
g = tsppi.getGraph()

clusterers = [ppi_networkit.PLP, ppi_networkit.PLM, ppi_networkit.PLM2, ppi_networkit.CNM]

clusterer = ppi_networkit.PLM(gamma=1000.0)
clusters = clusterer.run(g)

# TODO: do all this systematically, apparently lots of clusters are
# highly enriched
# THUS TODO: figure out how to continue from here and how to test multiple
#            clusterers
print("getting all genes")
all_genes = tsppi.getAllGenes()
print("loading go associations file")
assoc = load_go_associations(con, all_genes)
print("scoring clusters:")
score_clusters(clusters, tsppi, all_genes, assoc)


##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
