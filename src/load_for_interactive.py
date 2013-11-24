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


##############################
# enable tab completion
##############################
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
