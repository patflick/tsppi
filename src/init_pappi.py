'''
Configures all the data files path, and imports
all of them into the sqlite database.

@author: Patrick Flick
'''

import pappi
import os

# whether to use string-db or ccsb:
USE_STRINGDB=False

###############################
# Set the correct data folder:
###############################

# Note: 
# the database folder (DB_DATA_FOLDER) should always be an
# internal drive for performance purposes

# for easy of changing :)
AT_THE_LAB=True

#
if (AT_THE_LAB):
    # at the lab
    DATA_FOLDER = '/home/flick/dev/ppi/hpa/data/'
    CHANG_DATA_FOLDER = '/home/flick/dev/pappi/data/'
    DB_DATA_FOLDER = '/cygdrive/d/PPI/'
    STRINGDB_FILE='/cygdrive/d/PPI/string-db/human_protein.links.v9.0_700.csv'

else:
    # at home:
    DATA_FOLDER='/home/patrick/dev/bio/data/'
    CHANG_DATA_FOLDER = '/home/patrick/dev/bio/pappi/data/'
    DB_DATA_FOLDER = DATA_FOLDER
    # TODO string-db file, should also be on internal drive (large file)


##################################
# Set the path to all data files
##################################

HPA_FILE  = DATA_FOLDER + 'hpa_normal_tissue_v11.csv'
HGNC_FILE = DATA_FOLDER + 'hgnc_entrez_ensembl_uniprot.txt'
CCSB_FILE = DATA_FOLDER + 'HI_2012_PRE.tsv'
P2G_FILE = DATA_FOLDER + 'ensembl_ID_matching.csv'
BIOMART_FILE = DATA_FOLDER + 'mart_export.csv'
DATABASE  = DB_DATA_FOLDER + 'hpaDB.sqlite'

# load hk and ts genes
HK_FILE = CHANG_DATA_FOLDER + 'chang_hk.csv'
TS_FILE = CHANG_DATA_FOLDER + 'chang_ts.csv'

# output files
CCSB_PPI_OUT_FILE= DATA_FOLDER + 'ccsb_ppi.csv'
HPA_GENE_LEVELS_OUT_FILE=DATA_FOLDER + 'gene_levels.csv'



# TODO properly test this
# TODO continue with pipeline (merge two tissues, intersect with PPI, do graph scoring)
# TODO actually i can do graph scoring for _every_ tissue as pre-work
# and then comparing is as simple as loading and diffing two tables in the DB 


# first delete an old DB, to make sure everything is new
if (os.path.exists(DATABASE)):
    os.remove(DATABASE)

# get new database connection
con = pappi.sql.get_conn(DATABASE)

hpa_file = open(HPA_FILE)
print "Importing HPA data ..."
pappi.hpa.import_tissue(hpa_file, con)
print "Initializing HPA data ..."
pappi.hpa.init_gene_levels(con)

print "Importing Biomart table ..."
biomart_file = open(BIOMART_FILE)
pappi.matching.import_biomart_file(biomart_file, con)

if (USE_STRINGDB):
    ppi_file = open(STRINGDB_FILE)
    p2g_file = open(P2G_FILE)
    print "Importing string-db data ..."
    pappi.ppi.import_stringdb(ppi_file, p2g_file, con)
else:
    ccsb_file = open(CCSB_FILE)
    hgnc_file = open(HGNC_FILE)
    print "Importing CCSB data ..."
    pappi.ppi.import_ccsb(ccsb_file, hgnc_file, con)
    

# TODO: may need to load entrez to ensembl database first (in case string db is used)

print "Importing Chang HK/TS data ..."
hk_file = open(HK_FILE)
pappi.housekeeping.import_entrez_file(hk_file, con, "hk_entrez")
pappi.housekeeping.translate_entrez_2_ensembl(con, "hk_entrez", "hk_ensembl")

ts_file = open(TS_FILE)
pappi.housekeeping.import_entrez_file(ts_file, con, "ts_entrez")
pappi.housekeeping.translate_entrez_2_ensembl(con, "ts_entrez", "ts_ensembl")



ppi_out_file = open(CCSB_PPI_OUT_FILE, 'w')
print "Exporting PPI ..."
pappi.sql.dump_csv(ppi_out_file, "ppi_genes", con)

hpa_levels_file = open(HPA_GENE_LEVELS_OUT_FILE, 'w')
print "Exporting gene expression levels ..."
pappi.sql.dump_csv(hpa_levels_file, "hpa_gene_levels", con)

