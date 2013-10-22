import pappi.id_mapping
import pappi.sql
#import pappi.ppi_import
from pappi.ppi.ppi import PPI
from pappi.ppi.ccsb import CCSB
import os

DATA_FOLDER = '/home/patrick/dev/bio/data/'

HGNC_FILE = DATA_FOLDER + 'hgnc_entrez_ensembl_uniprot.txt'
BIOMART_FILE = DATA_FOLDER + 'mart_export.csv'
DATABASE = DATA_FOLDER + 'test_matching.sqlite'
CCSB_FILE = DATA_FOLDER + 'HI_2012_PRE.tsv'

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
#pappi.ppi_import.import_ccsb_file(ccsb_file, con, 'ccsb')

#pappi.id_mapping.map_identifier('ccsb', ['Gene_IDA', 'Gene_IDB'], 'entrez',
#                                'ccsb_hgnc1', 'hgnc', con, True)

ccsb_ppi = CCSB(CCSB_FILE, con)
ccsb_ppi.init_ppi(True)
