#!/usr/bin/env python3
#
# This script imports the data into the shared SQLite database
# The following data is imported and processed into unified formats:
# - ID mapping tables
# - GO terms associations
# - PPIs (5 major PPI networks):
#       * Dana Farber CCSB
#       * Bossi&Lehner
#       * StringDB
#       * PSICQUIC composite
#       * Havugimana et al.
# - Expression data sets:
#       * HPA
#       * Body Map 2.0
#       * RNAseq Atlas
#       * Gene Atlas
# - Expression `core` datasets
# - Gene + edge overlaps of datasets (PPI vs PPI, Expr vs Expr, PPI vs. Expr)
# - Expression aggregation statistics

import os
import re
import pappi.id_mapping
import pappi.sql
from pappi.data_config import *

##############################
# get database connection
##############################

# first delete an old DB, to make sure everything is new
if (os.path.exists(DATABASE)):
    os.remove(DATABASE)

# get new database connection
con = pappi.sql.get_conn(DATABASE)


##############################
# import the mapping tables
##############################

pappi.id_mapping.import_biomart_file(BIOMART_FILE, con)
pappi.id_mapping.import_hgnc_file(HGNC_FILE, con)


#####################
#  GO terms import  #
#####################

# load GO terms
import pappi.go_import

# import GO associations
pappi.go_import.import_go_association(GO_ASSOC_FILE, con)


##############################
# import PPIs
##############################

from pappi.ppis.ccsb import CCSB
from pappi.ppis.bossi_lehner import Bossi_Lehner
from pappi.ppis.havu import Havugimana
from pappi.ppis.string import StringDB
from pappi.ppis.psicquic import Psicquic
from pappi.ppis.psicquic_comb import PsicquicAll

ccsb_ppi = CCSB(CCSB_FILE, con)
ccsb_ppi.init_ppi(True)

bossi_ppi = Bossi_Lehner(BOSSI_FILE, con)
bossi_ppi.init_ppi(True)

havu_ppi = Havugimana(HAVU_FILE, con)
havu_ppi.init_ppi(True)

string_ppi = StringDB(STRING_FILE, con)
string_ppi.init_ppi(True)

# import all PSICQUIC networks
psicquic_ppis = []
psicquic_ppi_names = []
for pf in PSICQUIC_FILES:
    # extract the service name form the file name via regex
    service_name = re.match(r'.*_([a-zA-Z0-9-]+)\.tsv$', pf).group(1)
    # replace '-' (non SQL compatible) and set all to lower case
    service_name = service_name.replace('-', '_').lower()
    # create ppi class and init
    p = Psicquic(pf, con, service_name)
    p.init_ppi(True)
    # append to list of all ppis
    psicquic_ppis.append(p)
    psicquic_ppi_names.append(service_name)

# construct combined PSICQUIC PPI network
psicquic_all_ppi = PsicquicAll(con, psicquic_ppi_names)
psicquic_all_ppi.init_ppi(True)

##############################
# import expression data sets
##############################

from pappi.expr.hpa import HPA
from pappi.expr.hpa_all import HPA_All
from pappi.expr.emtab import Emtab
from pappi.expr.rnaseq_atlas import RnaSeqAtlas
from pappi.expr.gene_atlas import GeneAtlas

hpa_expr = HPA(HPA_FILE, con)
hpa_expr.init_data()

hpa_all_expr = HPA_All(HPA_FILE, con)
hpa_all_expr.init_data()

emtab_expr = Emtab(EMTAB_FILE, con)
emtab_expr.init_data()

rnaseq_atlas = RnaSeqAtlas(RNASEQ_ATLAS_FILE, con)
rnaseq_atlas.init_data()

gene_atlas = GeneAtlas(GENE_ATLAS_FILE, con)
gene_atlas.init_data()


##############################
# run overlap analyses
##############################

from pappi import overlap_analysis

overlap_analysis.calc_ppi_edge_overlap(con)
overlap_analysis.calc_ppi_id_overlap(con)
overlap_analysis.calc_expr_overlap(con)
# overlap of ppis and expression data sets (protein coverage of ppis by expr)
overlap_analysis.calc_pairwise_expr_ppi_id_overlap(con)
overlap_analysis.calc_pairwise_expr_ppi_edge_overlap(con)

# overlap of PPIs with each other (both IDs and edges)
overlap_analysis.calc_pairwise_ppi_id_overlap(con)
overlap_analysis.calc_pairwise_ppi_edge_overlap(con)


###################################
#  create core expression tables  #
###################################
for expr in [hpa_expr, hpa_all_expr, emtab_expr, rnaseq_atlas, gene_atlas]:
    print("creating core table for: " + expr.name)
    expr.create_core_table()
    # get expression counts for the core tables
    expr.expr_counts(True)
