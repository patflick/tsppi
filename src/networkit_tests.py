# this file is for interactive tests
# thus it just loads the configuration and
# database connection and enables tab completion

import os
import re
import itertools
import numpy

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

# import EnrichmentStudy
import sys
sys.path.append("/home/patrick/dev/bio/goatools")
#from goatools.go_enrichment import GOEnrichmentStudy


def enrichment_scoring(cluster, all_genes, assoc):
    print("got cluster of size " + str(len(cluster)) + " and #genes: "
          + str(len(all_genes)))
    print("Performing gene enrichment study:")
    population = set(all_genes)
    study = set(cluster)
    ge = GOEnrichmentStudy(population, assoc, obo_dag, alpha=0.05,
                           study=study,
                           methods=["bonferroni", "sidak", "holm"])
    ge.print_summary(pval=0.05)


def score_clusters(clusters, tsppi, sim_scorer):
    cluster_sizes = clusters.subsetSizeMap()
    # loop through all the clusters, get the associated gene names
    # and call the GO enrichment scoring of the cluster
    for cluster, size in cluster_sizes.items():
        cluster_nodes = clusters.getMembers(cluster)
        cl_node_names = [tsppi.getGeneName(x) for x in cluster_nodes]
        #enrichment_scoring(cl_node_names, all_genes, assoc)

        # calculate the cluster score
        in_score, ext_score = sim_scorer.gene_cluster_score(cl_node_names)
        # print it (TODO: do something with it (for statistics)
        print(str(cluster) + "\t" + str(size) + "\t" + str(in_score)
              + "\t" + str(ext_score))


#############################
#  get database connection  #
#############################

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

clusterers = [ppi_networkit.PLP, ppi_networkit.PLM, ppi_networkit.CNM]

clusterer = ppi_networkit.PLM(gamma=1000.0)
clusters = clusterer.run(g)

# import similarity scorer
from pappi.go.gene_prebuf_similarity import GoGenePreBufSimilarity

import time

start = time.time()
scorer = GoGenePreBufSimilarity(GO_OBO_FILE, GO_SCORE_FILE,
                                GO_SCORE_MAP_FILE, GO_BPSCORE_FILE,
                                GO_BPSCORE_MAP_FILE, con, True)
init_time = time.time() - start


print("scoring all clusters")
start = time.time()
score_clusters(clusters, tsppi, scorer)
score_time = time.time() - start

print("timing: init = " + str(init_time) + "s, score = "
      + str(score_time) + ", total = "
      + str(init_time + score_time))


##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
