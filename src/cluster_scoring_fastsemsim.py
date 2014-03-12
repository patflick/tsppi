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


##############################
# get database connection
##############################

con = pappi.sql.get_conn(DATABASE)


#######################################################
#  FastSemSim functions for GO functional similarity  #
#######################################################

gene1 = "EHF"
gene2 = "EZH2"
#gene2 = "HK1"

from pappi.go_similarity import GoSimilarity

# create go similarity scorer (based on fastSemSim)
go_sim = GoSimilarity(GO_OBO_FILE, GO_ASSOC_FILE, con, True)

print("BPscore:")
print(go_sim.pairwise_bp_score(gene1, gene2))

def cluster_get_n_go(cl_node_names):
    total_size = 0
    go_set = set()
    for gene in cl_node_names:
        if not gene in go_sim.gene_assoc:
            continue
        total_size = total_size + len(go_sim.gene_assoc[gene])
        go_set = go_set.union(go_sim.gene_assoc[gene])
    return (total_size, len(go_set))


def score_clusters(clusters, tsppi):
    cluster_sizes = clusters.subsetSizeMap()
    # loop through all the clusters, get the associated gene names
    # and call the GO enrichment scoring of the cluster
    print("scoring all clusters:")
    print("cl\tsize\tscore")
    for cluster, size in cluster_sizes.items():
        cluster_nodes = clusters.getMembers(cluster)
        cl_node_names = [tsppi.getGeneName(x) for x in cluster_nodes]
        #avg_score = go_sim.avg_bp_score(cl_node_names)
        avg_score = cluster_get_n_go(cl_node_names)
        print(str(cluster) + "\t" + str(size) + "\t" + str(avg_score))


#################################
#  Import PPI NetworKit module  #
#################################

#sys.path.append("/home/patrick/dev/bio/NetworKit-Flick/cython")
sys.path.append("/home/patrick/dev/bio/ppi_networkit/cython")
#import NetworKit
import ppi_networkit

# get new database connection
print("loading tsppi graph")

sqlio = ppi_networkit.SQLiteIO(DATABASE)

tsppi = sqlio.load_tsppi_graph("ccsb", "gene_atlas")
g = tsppi.getGraph()

clusterers = [ppi_networkit.PLP, ppi_networkit.PLM, ppi_networkit.PLM2, ppi_networkit.CNM]

#clusterer = ppi_networkit.PLM(gamma=1000.0)
clusterer = ppi_networkit.PLP()
clusters = clusterer.run(g)

# TODO: do all this systematically, apparently lots of clusters are
# highly enriched
# THUS TODO: figure out how to continue from here and how to test multiple
#            clusterers
print("scoring clusters:")
score_clusters(clusters, tsppi)


##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
