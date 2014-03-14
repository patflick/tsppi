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

import sys
sys.path.append("/home/patrick/dev/bio/goatools")
#from goatools.go_enrichment import GOEnrichmentStudy
from pappi.go_fastdag import GODag
from pappi.go_fastdag import name2id

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

def get_all_go_terms(assoc):
    terms = set()
    for ts in assoc.values():
        terms = terms.union(ts)
    return terms

def get_SimRel(go_dag, term1, term2):
    #print("terms: " + term1 + ", " + term2)
    term1 = name2id(term1)
    term2 = name2id(term2)
    LCAs = go_dag.get_lca(term1, term2)
    lca_IC = max(go_dag.IC[t] for t in LCAs)
    lca_p = min(go_dag.p[t] for t in LCAs)
    # TODO: check the formula (whether it is correct)
    denom = go_dag.IC[term1] + go_dag.IC[term2]
    if denom != 0:
        score = 2*lca_IC / denom * (1 - lca_p)
    else:
        score = 0.0
    if (score > 1.0 or score < 0.0):
        raise Exception
    return score

def get_bpscore(go_dag, assoc, gene1, gene2):
    terms1 = assoc[gene1]
    terms2 = assoc[gene2]
    m = -1
    score = 0.0
    for t1 in terms1:
        if not go_dag.has_term(t1):
            continue
        for t2 in terms2:
            if not go_dag.has_term(t2):
                continue
            score = get_SimRel(go_dag, t1, t2)
            score = max(m, score)
    if (score > 1.0 or score < 0.0):
        raise Exception
    return score


def avg_bp_score(go_dag, assoc, gene_set):
    """
    Calculates the average BPscore for a set of genes by first determining
    the BPscore between all pairs of genes in the set and then averaging
    over all BPscores.
    """
    genes = set()
    for g in gene_set:
        if g in assoc:
            genes.add(g)

    # get counts (merely for status print out)
    total_terms = sum(len(assoc[g]) for g in genes)
    unique_terms = set()
    for g in genes:
        unique_terms = unique_terms.union(assoc[g])
    print("scoring " + str(len(genes)) + " (" + str(len(gene_set)) + ") genes"
          " with " + str(total_terms) + " total and "
          + str(len(unique_terms)) + " unique terms")

    scores = []
    for g1, g2 in itertools.combinations(genes, 2):
        scores.append(get_bpscore(go_dag, assoc, g1, g2))
    # get mean
    avg = numpy.mean(scores)
    return avg

# load the GO Dag only once
obo_dag = GODag(obo_file=GO_OBO_FILE, only_namespace="biological_process")


def enrichment_scoring(cluster, all_genes, assoc):
    print("got cluster of size " + str(len(cluster)) + " and #genes: " + str(len(all_genes)))
    print("Performing gene enrichment study:")
    population = set(all_genes)
    study = set(cluster)
    ge = GOEnrichmentStudy(population, assoc, obo_dag, alpha=0.05,
            study=study, methods=["bonferroni", "sidak", "holm"])
    ge.print_summary(pval=0.05)


def score_clusters(clusters, tsppi, all_genes, assoc):
    cluster_sizes = clusters.subsetSizeMap()
    # TODO: do this globally
    print("get prob of nodes")
    obo_dag.term_probability(assoc)
    print("get IC")
    obo_dag.term_IC()
    # loop through all the clusters, get the associated gene names
    # and call the GO enrichment scoring of the cluster
    for cluster, size in cluster_sizes.items():
        cluster_nodes = clusters.getMembers(cluster)
        cl_node_names = [tsppi.getGeneName(x) for x in cluster_nodes]
        #enrichment_scoring(cl_node_names, all_genes, assoc)
        print("calc cluster of size: " + str(size))
        score = avg_bp_score(obo_dag, assoc, cl_node_names)
        print(str(cluster) + "\t" + str(size) + "\t" + str(score))


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

print("calc freq and probability")
obo_dag.term_probability(assoc)
print([obo_dag.freq[x] for x in obo_dag.roots])

print("scoring clusters:")
score_clusters(clusters, tsppi, all_genes, assoc)


##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
