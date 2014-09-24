# this file is for interactive tests of the ppi_networkit cython API
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


#############################
#  get database connection  #
#############################

# get new database connection
con = pappi.sql.get_conn(DATABASE)


#################################
#  Import PPI NetworKit module  #
#################################

sys.path.append(os.path.join(os.path.dirname(__file__), "../ppi_networkit/cython"))
#import NetworKit
import ppi_networkit
# set the log-level to "ERROR", this will ignore the [INFO] and [WARN] logs
ppi_networkit.setLogLevel("ERROR")
sqlio = ppi_networkit.SQLiteIO(DATABASE)

#ppi = "string"
#expr = "gene_atlas"
ppi = "ccsb"
ppi = "bossi"
expr = "hpa"
expr = "hpa_all"
tsppi = sqlio.load_tsppi_graph(ppi, expr)



##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
