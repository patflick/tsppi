import pappi.id_mapping
import pappi.sql
#import pappi.ppis_import
from pappi.ppis.ccsb import CCSB
from pappi.ppis.bossi_lehner import Bossi_Lehner
from pappi.ppis.havu import Havugimana
from pappi.expr.hpa import HPA
from pappi.expr.emtab import Emtab
import os

from pappi.data_config import *

# first delete an old DB, to make sure everything is new
if (os.path.exists(DATABASE)):
    os.remove(DATABASE)

# get new database connection
con = pappi.sql.get_conn(DATABASE)
pappi.id_mapping.import_biomart_file(BIOMART_FILE, con)
pappi.id_mapping.import_hgnc_file(HGNC_FILE, con)
#pappi.id_mapping.create_mapping_table("ensembl", "hgnc", con, True)
#pappi.id_mapping.create_mapping_table("uniprot", "hgnc", con, True)
#pappi.id_mapping.create_mapping_table("entrez", "hgnc", con, True)
#pappi.id_mapping.create_mapping_table("entrez", "uniprot", con, True)

# import CCSB to test the mapping

#ccsb_file = open(CCSB_FILE)
#print("Importing CCSB data ...")
#pappi.ppis_import.import_ccsb_file(ccsb_file, con, 'ccsb')

#pappi.id_mapping.map_identifier('ccsb', ['Gene_IDA', 'Gene_IDB'], 'entrez',
#                                'ccsb_hgnc1', 'hgnc', con, True)

ccsb_ppi = CCSB(CCSB_FILE, con)
ccsb_ppi.init_ppi(True)
#
#print("trying with bossi lehner PPI")
bossi_ppi = Bossi_Lehner(BOSSI_FILE, con)
##bossi_ppi.import_raw_file()
bossi_ppi.init_ppi(True)

#havu_ppi = Havugimana(HAVU_FILE, con)
#havu_ppi.init_ppi(True)
#
# try importing expression data sets
hpa_expr = HPA(HPA_FILE, con)
hpa_expr.init_data()

emtab_expr = Emtab(EMTAB_FILE, con)
emtab_expr.init_data()
