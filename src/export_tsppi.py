# this file is for interactive tests
# thus it just loads the configuration and
# database connection and enables tab completion

import os
import re
import pappi.id_mapping
import pappi.sql
from pappi.data_config import *

##############################
# get database connection
##############################

# get new database connection
con = pappi.sql.get_conn(DATABASE)


# for testing expr and ppis (bossi and hpa as example)
from pappi.expr.hpa import *
from pappi.ppis.bossi_lehner import *
# create classes (don't init data)
hpa = HPA("", con)
bossi = Bossi_Lehner("", con)
# get simple variables
e = hpa
p = bossi

NODE_LABELS_FILE = '/home/patrick/dev/bio/test_data/node_labels.txt'
EDGE_LIST_FILE = '/home/patrick/dev/bio/test_data/edge_list.txt'

# create export tables, that are read by the C++ SQLite part for creating of
# the tissue specific graph for NetworKit
e.export_node_labels(p.name + '_ids')
p.export_to_edge_list(only_ids=e.name + '_ids')


##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
