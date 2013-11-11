import os
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


##############################
# import PPIs
##############################

from pappi.ppis.ccsb import CCSB
from pappi.ppis.bossi_lehner import Bossi_Lehner
from pappi.ppis.havu import Havugimana

ccsb_ppi = CCSB(CCSB_FILE, con)
ccsb_ppi.init_ppi(True)

bossi_ppi = Bossi_Lehner(BOSSI_FILE, con)
bossi_ppi.init_ppi(True)

havu_ppi = Havugimana(HAVU_FILE, con)
havu_ppi.init_ppi(True)


##############################
# import expression data sets
##############################

from pappi.expr.hpa import HPA
from pappi.expr.emtab import Emtab
from pappi.expr.rnaseq_atlas import RnaSeqAtlas
from pappi.expr.gene_atlas import GeneAtlas

hpa_expr = HPA(HPA_FILE, con)
hpa_expr.init_data()

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

# phyper(8019, n*(n-1)/2 - k*(k-1)/2, k*(k-1)/2,13943)
