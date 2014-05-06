# this file is for interactive tests
# thus it just loads the configuration and
# database connection and enables tab completion

import os
import re
import itertools
import numpy
import time

import pappi.id_mapping
import pappi.sql
from pappi.data_config import *

##############################################################
#  TODO: temporary content (needs to be put into a module)   #
##############################################################

# histogram function using counter
from collections import Counter

# import EnrichmentStudy
import sys
sys.path.append("/home/patrick/dev/bio/goatools")
#from goatools.go_enrichment import GOEnrichmentStudy
# import PPI networkit
#sys.path.append("/home/patrick/dev/bio/NetworKit-Flick/cython")
sys.path.append("/home/patrick/dev/bio/ppi_networkit/cython")
#import NetworKit
import ppi_networkit
# set the log-level to "ERROR", this will ignore the [INFO] and [WARN] logs
ppi_networkit.setLogLevel("ERROR")


# TODO: put these somewhere unified
PPI_NAMES = ["string", "ccsb", "bossi", "psicquic_all", "havu"]
EXPR_NAMES = ["hpa", "hpa_all", "emtab", "rnaseq_atlas", "gene_atlas"]


def cluster_hist(clusters):
    c = Counter()
    cluster_sizes = clusters.subsetSizes()
    for size in cluster_sizes:
        c[size] += 1
    return c


def agg_cluster_size_hist(clusters, bins=5):
    cluster_sizes = clusters.subsetSizes()
    max_size = max(cluster_sizes)
    hist = [0]*bins
    # bin size as ceiling of max_size / bin
    bin_size = (max_size-1) // bins + 1
    for size in cluster_sizes:
        b = (size-1) // bin_size
        hist[b] += 1

    # output accumulated histogram of sizes
    print("\t| ".join("<" + str(i * bin_size) for i in range(1, bins+1)))
    print("\t| ".join(str(n_bin) for n_bin in hist))


def clusters_weighted_mean(clusters):
    cluster_sizes = clusters.subsetSizes()
    n = sum(cluster_sizes)
    # the size is also the mean, therefore sum over size^2
    weighted_mean = sum(s*s for s in cluster_sizes) / n
    return weighted_mean


def clusters_weighted_median(clusters):
    cluster_sizes = clusters.subsetSizes()
    n = sum(cluster_sizes)
    med = n / 2
    cluster_sizes = sorted(cluster_sizes)
    s = 0
    for size in cluster_sizes:
        s += size
        if s >= med:
            return size
    return max(cluster_sizes)


def print_sorted_hist(hist):
    print("size\t#clusters")
    for size, num in iter(sorted(hist.items())):
        print(str(size) + ":\t" + str(num))



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
    scores = []
    for cluster, size in cluster_sizes.items():
        cluster_nodes = clusters.getMembers(cluster)
        cl_node_names = [tsppi.getGeneName(x) for x in cluster_nodes]
        #enrichment_scoring(cl_node_names, all_genes, assoc)

        # calculate the cluster score
        if len(cl_node_names) >= 2:
            in_score, ext_score = sim_scorer.gene_cluster_score(cl_node_names)
            scores.append((cluster, size, in_score - ext_score))
        else:
            scores.append((cluster, size, 0.0))

    return scores


def clusters_modul(clusters, graph):
    modsPerCluster = ppi_networkit.ModularityPerCluster()
    mods = modsPerCluster.getModularities(clusters, graph)
    return mods


#############################
#  get database connection  #
#############################

# get new database connection
con = pappi.sql.get_conn(DATABASE)


class StdoutWriter:
    def __init__(self):
        pass

    def writerow(self, row):
        print(row)


class SQLWriter:
    def __init__(self, con, drop_table=True):
        self.con = con
        self.ppi = ""
        self.expr = ""
        self.cur = con.cursor()
        if drop_table:
            self.cur.execute('DROP TABLE IF EXISTS clustering_scoring_results')
        self.cur.execute('CREATE TABLE IF NOT EXISTS clustering_scoring_results '
                         '(ppi, expr, clusterer, type, cluster_id, size, '
                         ' bpscore, modularity)')

    def set_ppi(self, ppi):
        self.ppi = ppi

    def set_expr(self, ppi):
        self.expr = expr

    def set_clusterer(self, clusterer):
        self.clusterer = clusterer

    def writerow(self, row):
        self.cur.execute('INSERT INTO clustering_scoring_results '
                         'VALUES (?,?,?,?,?,?,?,?)',
                         tuple([self.ppi, self.expr, self.clusterer]
                               + list(row)))

    def commit(self):
        """"
        Commits all current changes and re-opens the cursor.
        """
        self.cur.close()
        self.con.commit()
        self.cur = self.con.cursor()

    def __del__(self):
        self.cur.close()
        self.con.commit()


#################################
#  Import PPI NetworKit module  #
#################################

def run_and_score_clustering(graph, clusterer, scorer, writer=None, category=None):
    # run clustering
    start = time.time()
    clusters = clusterer.run(g)
    t_cluster = time.time() - start
    # print histogram
    #hist = cluster_hist(clusters)
    #print_sorted_hist(hist)
    #agg_cluster_size_hist(clusters)
    weighted_mean = clusters_weighted_mean(clusters)
    w_median = clusters_weighted_median(clusters)

    # run scoring
    start = time.time()
    scores = score_clusters(clusters, tsppi, scorer)
    t_score = time.time() - start

    # cluster modularities
    mods = clusters_modul(clusters, graph)

    avg_score = numpy.mean([score for i, size, score in scores])

    # write out result
    if not writer is None:
        for i, size, score in scores:
            writer.writerow([category, i, size, score, mods[i]])

    # print some scoring results:
    print("avg score: " + str(avg_score) + ", time clustering: "
          + str(t_cluster) + ", time scoring: " + str(t_score)
          + " wmean: " + str(weighted_mean)
          + " wmedian: " + str(w_median))


def run_multiple_clusterers(graph, scorer):
    clusterers = [ppi_networkit.PLP, ppi_networkit.PLM, ppi_networkit.CNM]
    for gamma in [1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 100.0, 200.0, 500.0, 1000.0]:
        clusterer = ppi_networkit.PLM(gamma=gamma)
        print("PLM with gamma = " + str(gamma))
        run_and_score_clustering(graph, clusterer, scorer)
    print("===================================================")
    clusterer = ppi_networkit.PLP()
    run_and_score_clustering(graph, clusterer, scorer)


def run_ts_clustering(tsppi, clusterer, scorer, writer=None):
    print("GLOBAL GRAPH:")
    g = tsppi.getGraph()
    run_and_score_clustering(g, clusterer, scorer, writer, "Global")

    print("TS GRAPHs")
    nTissues = tsppi.getNumberOfTissues()

    for t in range(0, nTissues):
        tissue_name = tsppi.getTissueName(t)
        print("Tissue: " + tissue_name)
        ts_graph = tsppi.getTsGraph(t)
        run_and_score_clustering(ts_graph, clusterer, scorer, writer, tissue_name)


def run_edgescore_clustering(tsppi, clusterer, scorer, writer=None):
    print("Scoring on EdgeScore graph (TS/Global hybrid via edge weighting)")
    g = tsppi.getEdgeScoreGraph()
    run_and_score_clustering(g, clusterer, scorer, writer, "EdgeScoring")


def run_global_clustering(tsppi, clusterer, scorer, writer=None):
    print("scoring on global graph:")
    g = tsppi.getGraph()
    run_and_score_clustering(g, clusterer, scorer, writer, "GLOBAL")


# import similarity scorer
from pappi.go.gene_prebuf_similarity import GoGenePreBufSimilarity


def get_scorer(con):
    start = time.time()
    scorer = GoGenePreBufSimilarity(GO_OBO_FILE, GO_SCORE_FILE,
                                    GO_SCORE_MAP_FILE, GO_BPSCORE_FILE,
                                    GO_BPSCORE_ROW_FILE, GO_BPSCORE_MAP_FILE,
                                    con, True)
    init_time = time.time() - start
    print("scorer init time: " + str(init_time) + " s")
    return scorer


sqlio = ppi_networkit.SQLiteIO(DATABASE)

# initialize GO/BP scorer
scorer = get_scorer(con)

# create results writer
writer = SQLWriter(con)

# get clusterer
#clusterers = [ppi_networkit.PLP, ppi_networkit.PLM, ppi_networkit.CNM]
gamma = 1
clusterer = ppi_networkit.PLM(gamma=gamma)
writer.set_clusterer("PLM-gamma-1.0")

for ppi in PPI_NAMES:
    for expr in EXPR_NAMES:
        print()
        print("##################################################")
        print(" getting graph properties of `" + ppi + "_" + expr + "`")
        print("##################################################")

        # set the current ppi and expr
        writer.set_ppi(ppi)
        writer.set_expr(expr)

        # get graph
        tsppi = sqlio.load_tsppi_graph(ppi, expr)
        g = tsppi.getGraph()

        # run the clustering algos
        run_global_clustering(tsppi, clusterer, scorer, writer)
        run_ts_clustering(tsppi, clusterer, scorer, writer)
        run_edgescore_clustering(tsppi, clusterer, scorer, writer)

        # commit all current changes to the SQL server
        writer.commit()


##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
