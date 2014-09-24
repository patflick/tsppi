#!/usr/bin/env python3
#
# Uses clustering/community detection algorithms to find functionally related
# modules in tissue-specific PPI networks.
# - executes PLM, PLP and CNM clusterers from NetworKit on different graphs
# - saves modularities and per cluster modularities and BPScores into SQL table
#   for later analysis

import os
import re
import itertools
import numpy
import time

import pappi.id_mapping
import pappi.sql
from pappi.data_config import *

# histogram function using counter
from collections import Counter


# import PPI networkit
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../ppi_networkit/cython"))
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


def score_clusters(clusters, tsppi, sim_scorer):
    cluster_sizes = clusters.subsetSizeMap()
    # loop through all the clusters, get the associated gene names
    # and call the GO enrichment scoring of the cluster
    scores = []
    for cluster, size in cluster_sizes.items():
        cluster_nodes = clusters.getMembers(cluster)
        cl_node_names = [tsppi.getGeneName(x) for x in cluster_nodes]

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
        self.cur.execute('CREATE TABLE IF NOT EXISTS clustering_scoring_results'
                         ' (ppi, expr, clusterer, type, cluster_id, size, '
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

def run_and_score_clustering(graph, clusterer, scorer,
                             writer=None, category=None):
    # run clustering
    start = time.time()
    clusters = clusterer.run(graph)
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
    for gamma in [1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0,
                  50.0, 100.0, 200.0, 500.0, 1000.0]:
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
        run_and_score_clustering(ts_graph, clusterer,
                                 scorer, writer, tissue_name)


def run_edgescore_clustering(tsppi, clusterer, scorer, writer=None):
    print("Scoring on EdgeScore graph (TS/Global hybrid via edge weighting)")
    #g = tsppi.getEdgeCorrelationGraph()
    #run_and_score_clustering(g, clusterer, scorer, writer, "EdgeCorrelation")
    g = tsppi.getEdgeCoexprCountGraph()
    run_and_score_clustering(g, clusterer, scorer, writer, "EdgeCoexprCount")


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


if __name__ == "__main__":

    sqlio = ppi_networkit.SQLiteIO(DATABASE)

    # initialize GO/BP scorer
    scorer = get_scorer(con)

    # create results writer
    writer = SQLWriter(con, True)

    # get clusterer
    #clusterers = [ppi_networkit.PLP, ppi_networkit.PLM, ppi_networkit.CNM]
    clusterers = [("PLM-gamma-5.0", ppi_networkit.PLM(gamma=5)),
                  ("PLM-gamma-10.0", ppi_networkit.PLM(gamma=10)),
                  ("PLM-gamma-50.0", ppi_networkit.PLM(gamma=50)),
                  ("PLP", ppi_networkit.PLP()),
                  ("CNM", ppi_networkit.CNM())]
    #clusterers = [("PLM-gamma-100.0", ppi_networkit.PLM(gamma=100))]
    #clusterers = [("PLM-gamma-50.0", ppi_networkit.PLM(gamma=50))]

    # run each clusterer on all ppis and expression datasets
    for clusterer_name, clusterer in clusterers:
        writer.set_clusterer(clusterer_name)

        for ppi in PPI_NAMES:
            for expr in EXPR_NAMES:
                print()
                print("##################################################")
                print(" getting graph properties of " + ppi + "_" + expr)
                print("##################################################")

                # set the current ppi and expr
                writer.set_ppi(ppi)
                writer.set_expr(expr)

                # get graph
                tsppi = sqlio.load_tsppi_graph(ppi, expr)

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
