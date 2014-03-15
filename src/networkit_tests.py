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

def all_term_sim(go_dag, assoc):
    # get the unique go_terms
    terms = set()
    for term_set in assoc.values():
        terms = terms.union(term_set)
    # intersect with terms in GO Dag
    terms = terms.intersection(go_dag.terms)
    
    # number of terms
    nTerms = len(terms)

    # map terms to numberical range [0, nTerms-1]
    i = 0
    term_mapping = dict()
    idx_2_term = list()
    for t in terms:
        term_mapping[t] = i
        idx_2_term.append(t)
        i = i + 1

    print("Pre-calculating SemSim between all " + str(len(terms)) + " terms.")
    sim_vals = numpy.zeros((nTerms, nTerms))
    
    # fill matrix in row major
    for i in range(0, nTerms):
        t1 = idx_2_term[i]
        if i % 10 == 0:
            print("filling row " + str(i) + "/" + str(nTerms))
        for j in range(i+1, nTerms):
            t2 = idx_2_term[j]
            sim = get_SimRel(go_dag, t1, t2)
            sim_vals[i][j] = sim
            sim_vals[j][i] = sim
    numpy.save("half-simrel", sim_vals)

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
        print("calc cluster of size: " + str(size))
        score = sim_scorer.gene_set_score(cl_node_names)
        #score = avg_bp_score(obo_dag, assoc, cl_node_names)
        print(str(cluster) + "\t" + str(size) + "\t" + str(score))


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

clusterers = [ppi_networkit.PLP, ppi_networkit.PLM, ppi_networkit.PLM2,
              ppi_networkit.CNM]

clusterer = ppi_networkit.PLM(gamma=1000.0)
#clusters = clusterer.run(g)

# TODO: do all this systematically, apparently lots of clusters are
# highly enriched
# THUS TODO: figure out how to continue from here and how to test multiple
#            clusterers
print("getting all genes")
all_genes = tsppi.getAllGenes()

from pappi.go_fast_similarity import GoFastSimilarity
from pappi.go_fastSemSim_similarity import GoFastSemSimSimilarity

# import my GO scoring
simScorer = GoFastSimilarity(GO_OBO_FILE, con, True)


#all_term_sim(obo_dag, assoc)

gene1 = "EHF"
gene2 = "EZH2"
#gene2 = "HK3"
print("BPscore:")
print(simScorer.gene_pairwise_score(gene1, gene2))

print("scoring clusters:")
#score_clusters(clusters, tsppi, all_genes, assoc)


##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
