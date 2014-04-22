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


###################################
#  get and save graph properties  #
###################################

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
    #save_node_properties(ppi_name, ppi.getAllGenes(), props, con)
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


def get_graph_properties(ppi, con, verbose=True, timings={}):
    # get graph and name
    g = ppi.getGraph()
    ppi_name = ppi.getPpiName()

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
    #save_graph_properties(ppi_name, props, con)
    return timings


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
for ppi_name in ["string", "ccsb", "bossi", "psicquic_all", "havu"]:
    print()
    print("##################################################")
    print(" getting graph properties of `" + ppi_name + "`")
    print("##################################################")
    ppi = sqlio.load_ppi_graph(ppi_name)
    timings = get_node_properties(ppi, con, True, timings)
    timings = get_graph_properties(ppi, con, True, timings)

print()
print("##################################################")
print("         time taken by each part")
print("##################################################")
print(timings)

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
