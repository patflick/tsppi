# this file is for interactive tests
# thus it just loads the configuration and
# database connection and enables tab completion

import os
import re
import itertools
import numpy
import time
import sys

import pappi.id_mapping
import pappi.sql
from pappi.data_config import *

# whether to save the calculated results
save_results = True

# TODO: do this kind of analysis later (nodes with higher betweenness as before)
#SELECT * FROM ccsb_hpa_ts_node_properties AS a INNER JOIN hpa as B ON a.Gene =
# b.Gene AND a.Tissue = b.Type INNER JOIN ccsb_node_properties AS c ON a.Gene =
# c.Gene WHERE b.Expressed = 1 AND a.Betweenness > c.Betweenness

# TODO: put these somewhere unified
PPI_NAMES = ["string", "ccsb", "bossi", "psicquic_all", "havu"]
EXPR_NAMES = ["hpa", "hpa_all", "emtab", "rnaseq_atlas", "gene_atlas"]

# FIXME: limited testing set (that doesn't need 12h)
#PPI_NAMES = ["ccsb", "havu"]
#EXPR_NAMES = ["hpa", "rnaseq_atlas"]

###################################
#  get and save graph properties  #
###################################

# TODO: abstract and do the get/save for node and graph properties for all of
#       these classes of graphs:
#       - global graph (ppi, property)
#       - global ts graph (ppi, expr, property)
#       - all the ts graphs (ppi, expr, tissue, property) [actually put each ppi X expr into own table !?]


# save the properties from the TS graphs (PxExT)
def save_ts_node_properties(graph_name, gene_names, ts_names, properties, sql_conn):
    # initialize the cursor object
    cur = sql_conn.cursor()

    # create the properties table
    table = graph_name +  "_ts_node_properties"
    props = list(properties.keys())
    cur.execute('DROP TABLE IF EXISTS "' + table + '"')
    cur.execute('CREATE TABLE "' + table + '" ("Gene" varchar(16), '
                + '"Tissue" varchar(16), '
                + ', '.join('"' + p + '" REAL' for p in props) + ');')

    # fill the properties table
    for t in range(0, len(ts_names)):
        for g in range(0, len(gene_names)):
            prop_fields = ', '.join(props)
            values = tuple([gene_names[g], ts_names[t]] + [properties[p][t][g] for p in props])
            cur.execute('INSERT INTO "' + table + '" (Gene, Tissue, '
                        + prop_fields
                        + ') VALUES (?,?,' + ','.join('?'*len(props)) + ')',
                        values)

    # cluse cursor and commit changes
    cur.close()
    sql_conn.commit()

def get_ts_node_properties(tsppi, con, verbose=True, timings={}):
    # get name
    graph_name = tsppi.getPpiName()

    props = dict()
    if verbose:
        print("getting degrees")
    t = time.time()
    props['degree'] = tsppi.getTsDegrees()
    timings[(ppi_name,'ts_degrees')] = time.time() - t

    g = tsppi.getGraph()

    # calculate betweenness
    if verbose:
        print("getting graph betweenness")
    t = time.time()
    betweenness_scores = tsppi.getTsBetweenness()
    props['Betweenness'] = betweenness_scores
    timings[(ppi_name,'ts_betweenness')] = time.time() - t

    # calculate clustering coeff
    if verbose:
        print("getting graph clustering coeff")
    # TODO: clustering coeff for Ts graphs

    # save the node properties to the database
    if verbose:
        print("saving node properties")
    if save_results:
        save_ts_node_properties(graph_name, tsppi.getAllGenes(), tsppi.getAllTissues(), props, con)
    return timings


def save_node_properties(graph_name, gene_names, properties, sql_conn):
    for prop_name, props in properties.items():
        if len(gene_names) != len(props):
            raise Exception()

    # initialize the cursor object
    cur = sql_conn.cursor()

    # create the properties table
    table = graph_name + "_node_properties"
    props = list(properties.keys())
    cur.execute('DROP TABLE IF EXISTS "' + table + '"')
    cur.execute('CREATE TABLE "' + table + '" ("Gene" varchar(16), '
                + ', '.join('"' + p + '" REAL' for p in props) + ');')

    # fill the properties table
    for i in range(0, len(gene_names)):
        prop_fields = ', '.join(props)
        values = tuple([gene_names[i]] + [properties[p][i] for p in props])
        cur.execute('INSERT INTO "' + table + '" (Gene, ' + prop_fields
                    + ') VALUES (?,' + ','.join('?'*len(props)) + ')',
                    values)

    # cluse cursor and commit changes
    cur.close()
    sql_conn.commit()


def get_node_properties(ppi, con, verbose=True, timings={}):
    # get name
    ppi_name = ppi.getPpiName()

    props = dict()
    if verbose:
        print("getting degrees")
    t = time.time()
    props['degree'] = ppi.getDegrees()
    timings[(ppi_name,'degrees')] = time.time() - t

    g = ppi.getGraph()

    # calculate betweenness
    if verbose:
        print("getting graph betweenness")
    t = time.time()
    between = ppi_networkit.Betweenness(g)
    between.run()
    props['Betweenness'] = between.scores()
    timings[(ppi_name,'betweenness')] = time.time() - t

    # calculate clustering coeff
    if verbose:
        print("getting graph clustering coeff")
    t = time.time()
    clustercoef = ppi_networkit.ClusteringCoefficient()
    props['ClusteringCoeff'] = clustercoef.exactLocal(g)
    timings[(ppi_name, 'clusteringcoeff')] = time.time() - t

    # save the node properties to the database
    if verbose:
        print("saving node properties")
    if save_results:
        save_node_properties(ppi_name, ppi.getAllGenes(), props, con)
    return timings


def save_graph_properties(graph_name, props, sql_conn):
    # initialize the cursor object
    cur = sql_conn.cursor()

    # create the properties table
    table = graph_name + "_properties"
    cur.execute('DROP TABLE IF EXISTS "' + table + '"')
    cur.execute('CREATE TABLE "' + table + '" ("Property" varchar(16), '
                '"Value" REAL)')

    # for every property:
    for p, val in props.items():
        cur.execute('INSERT INTO "' + table + '" (Property,Value) '
                    'VALUES (?,?)', (p, val))
    # close connection and commit changes
    cur.close()
    sql_conn.commit()

def graph_stats(g, verbose=True, timings={}):
    props = dict()

    if verbose:
        print("get basic graph info")
    t = time.time()
    props['n'] = g.numberOfNodes()
    props['m'] = g.numberOfEdges()
    timings[(ppi_name, "size")] = time.time() - t

    if verbose:
        print("get degree stats")
    t = time.time()
    gp = ppi_networkit.GraphProperties
    minmax = gp.minMaxDegree(g)
    props['min_deg'] = minmax[0]
    props['max_deg'] = minmax[1]
    props['avg_deg'] = gp.averageDegree(g)
    props['deg_assortativity'] = gp.degreeAssortativity(g, False)
    timings[(ppi_name, "degree")] = time.time() - t

    if verbose:
        print("get connected components")
    t = time.time()
    cc = ppi_networkit.ConnectedComponents()
    cc.run(g)
    props['conn_comp'] = cc.numberOfComponents()
    timings[(ppi_name, "conn_comp")] = time.time() - t

    if verbose:
        print("get global clustering coeff")
    t = time.time()
    clc = ppi_networkit.ClusteringCoefficient()
    props['gl_cluster_coeff'] = clc.exactGlobal(g)
    timings[(ppi_name, "gl_cluster_coeff")] = time.time() - t


    if verbose:
        print("get diameter")
    t = time.time()
    if props['conn_comp'] == 1:
        dm = ppi_networkit.Diameter
        props['diameter'] = dm.exactDiameter(g)
    else:
        print("    ! graph not connected -> setting diameter = 0")
        props['diameter'] = 0
    timings[(ppi_name, "diameter")] = time.time() - t

    return (props, timings)

# TODO: this for tsppi subgraphs
def get_graph_properties(ppi, con, verbose=True, timings={}):
    # get graph and name
    g = ppi.getGraph()
    ppi_name = ppi.getPpiName()

    # get the graph statistics
    (props, timings) = graph_stats(g, verbose, timings)

    if save_results:
        save_graph_properties(ppi_name, props, con)
    return timings

def prepare_ts_graph_property_table(sql_conn):
    # initialize the cursor object
    cur = sql_conn.cursor()
    # create the properties table
    table = "ts_graph_properties"
    cur.execute('DROP TABLE IF EXISTS "' + table + '"')
    cur.execute('CREATE TABLE "' + table + '" ("ppi" varchar(16), "expr" varchar(16), "Tissue" varchar(32), "Property" varchar(16), '
                '"Value" REAL)')
    # close connection and commit changes
    cur.close()
    sql_conn.commit()

def save_ts_graph_properties(ppi_name, expr_name, tissue_name, props, sql_conn):
    # initialize the cursor object
    cur = sql_conn.cursor()
    table = "ts_graph_properties"

    # for every property:
    for p, val in props.items():
        cur.execute('INSERT INTO "' + table + '" (ppi,expr,Tissue,Property,Value) '
                    'VALUES (?,?,?,?,?)', (ppi_name, expr_name, tissue_name, p, val))
    # close connection and commit changes
    cur.close()
    sql_conn.commit()

# TODO: save global ts properties into shared table (ppi, expr, propr, value)
def get_ts_graph_properties(tsppi, ppi_name, expr_name, con, verbose=True,
                            timings={}):
    # for all tissues, get the tissue graph and then get and save the properties
    for t in range(0, tsppi.getNumberOfTissues()):
        # get tissue name
        tissue_name = tsppi.getTissueName(t)
        # get the tissue specific graph for the current tissue
        g = tsppi.getTsGraph(t)

        # get global properties of that graph
        (props, timings) = graph_stats(g, verbose, timings)

        if save_results:
            save_ts_graph_properties(ppi_name, expr_name, tissue_name, props, con)


# save all the timings into the SQL database
def save_timings(timings, sql_conn, table):
    # create new, empty table
    cur = sql_conn.cursor()
    cur.execute('DROP TABLE IF EXISTS "' + table + '"')
    cur.execute('CREATE TABLE "' + table + '" ("graph_name" varchar(16), "property" varchar(16), "time" REAL)')

    # insert all values
    for k,v in timings.items():
        (graph_name, prop) = k
        cur.execute('INSERT INTO "' + table + '" (graph_name,property,time) '
                    'VALUES (?,?,?)', (graph_name, prop, v))
    # close connection and commit changes
    cur.close()
    sql_conn.commit()



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
# set the log-level to "ERROR", this will ignore the [INFO] and [WARN] logs
ppi_networkit.setLogLevel("ERROR")
sqlio = ppi_networkit.SQLiteIO(DATABASE)


#################################
#  Get global-graph properties  #
#################################

timings = {}
for ppi_name in PPI_NAMES:
    print()
    print("##################################################")
    print(" getting graph properties of `" + ppi_name + "`")
    print("##################################################")
    ppi = sqlio.load_ppi_graph(ppi_name)
    timings = get_graph_properties(ppi, con, True, timings)
    timings = get_node_properties(ppi, con, True, timings)
    save_timings(timings, con, "ppi_graph_stats_timings")

##########################################
#  Get tissue-specific graph properties  #
##########################################

timings = {}
for ppi_name in PPI_NAMES:
    for expr_name in EXPR_NAMES:
        graph_name = ppi_name + "_" + expr_name
        print()
        print("##################################################")
        print(" getting graph properties of `" + graph_name + "`")
        print("##################################################")
        tsppi = sqlio.load_tsppi_graph(ppi_name, expr_name)

        # global ts graph properties
        timings = get_graph_properties(tsppi, con, True, timings)
        timings = get_node_properties(tsppi, con, True, timings)
        save_timings(timings, con, "ppi_graph_stats_timings")

sys.exit(1)

###########################################
#  Get tissue specific tissue properties  #
###########################################

timings = {}
prepare_ts_graph_property_table(con)
for ppi_name in PPI_NAMES:
    for expr_name in EXPR_NAMES:
        graph_name = ppi_name + "_" + expr_name
        print()
        print("##################################################")
        print(" getting TS tissue properties of `" + graph_name + "`")
        print("##################################################")
        tsppi = sqlio.load_tsppi_graph(ppi_name, expr_name)

        # get tissue specific graph properties (no timing for this available)
        get_ts_graph_properties(tsppi, ppi_name, expr_name, con, True)
        # get tissue specific node properties
        timings = get_ts_node_properties(tsppi, con, True, timings)


##########################################
#  Get tissue-specific graph properties  #
##########################################

#ppi = "ccsb"
#ppi = "bossi"
#expr = "hpa"
#expr = "hpa_all"
#tsppi = sqlio.load_tsppi_graph(ppi, expr)
#g = tsppi.getGraph()
#

################################
#  enable for interactive use  #
################################

import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
