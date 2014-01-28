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

# import stuff
import pappi.overlap_analysis

import pappi.expr
from pappi.ppis.bossi_lehner import Bossi_Lehner
# create classes (don't init data)
# for testing expr and ppis (bossi and hpa as example)
hpa = pappi.expr.HPA("", con)
bossi = Bossi_Lehner("", con)
# get simple variables
e = hpa
p = bossi

expr_classes = [pappi.expr.HPA, pappi.expr.HPA_All, pappi.expr.RnaSeqAtlas, pappi.expr.GeneAtlas, pappi.expr.Emtab]
exprs = [C("", con) for C in expr_classes]

emtab = pappi.expr.Emtab(EMTAB_FILE, con)

##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
