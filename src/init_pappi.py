'''
Configures all the data files path, and imports
all of them into the sqlite database.

@author: Patrick Flick
'''

import pappi
import os

# whether to use string-db or ccsb:
#PPI_TO_USE="string-db"
USE_STRINGDB_PPI=False
USE_CCSB_PPI=True
USE_MMC_PPI=False


DUMP_PPI=False
DUMP_GENE_EXPR_LEVELS=False

###############################
# Set the correct data folder:
###############################

# Note: 
# the database folder (DB_DATA_FOLDER) should always be an
# internal drive for performance purposes

# for easy of changing :)
AT_THE_LAB=False

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
MMC_FILE = CHANG_DATA_FOLDER + 'cell_havugimana_ppi.tsv'
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


print "Importing and initializing ID matching/mapping ..."
hgnc_file = open(HGNC_FILE)
biomart_file = open(BIOMART_FILE)
pappi.matching.import_mappings(hgnc_file, biomart_file, con)


print "Importing HPA data ..."
hpa_file = open(HPA_FILE)
pappi.hpa.import_tissue(hpa_file, con)
print "Initializing HPA data ..."
pappi.hpa.init_gene_levels(con)




if (USE_STRINGDB_PPI):
    ppi_file = open(STRINGDB_FILE)
    p2g_file = open(P2G_FILE)
    print "Importing string-db data ..."
    pappi.ppi.import_stringdb(ppi_file, p2g_file, con)

if (USE_CCSB_PPI):
    ccsb_file = open(CCSB_FILE)
    print "Importing CCSB data ..."
    pappi.ppi.import_ccsb(ccsb_file, hgnc_file, con)

if (USE_MMC_PPI):
    mmc_file = open(MMC_FILE)
    print "Importing MMC PPI data ..."
    pappi.ppi.import_mmc(mmc_file, con)
    

# after importing, initializing and transforming PPIs to use HGNC symbols
# continue with processing:
pappi.ppi.init_edge_expression(con)

# TODO: may need to load entrez to ensembl database first (in case string db is used)

if (False):
    print "Importing Chang HK/TS data ..."
    hgnc_file = open(HGNC_FILE)
    pappi.matching.import_hgnc_entrez2ensembl(hgnc_file, con)
    
    hk_file = open(HK_FILE)
    pappi.housekeeping.import_entrez_file(hk_file, con, "hk_entrez")
    pappi.housekeeping.translate_entrez_2_ensembl(con, "hk_entrez", "hk_ensembl")
    
    ts_file = open(TS_FILE)
    pappi.housekeeping.import_entrez_file(ts_file, con, "ts_entrez")
    pappi.housekeeping.translate_entrez_2_ensembl(con, "ts_entrez", "ts_ensembl")


if (DUMP_PPI):
    ppi_out_file = open(CCSB_PPI_OUT_FILE, 'w')
    print "Exporting PPI ..."
    pappi.sql.dump_csv(ppi_out_file, "ppi_genes", con)

if (DUMP_GENE_EXPR_LEVELS):
    hpa_levels_file = open(HPA_GENE_LEVELS_OUT_FILE, 'w')
    print "Exporting gene expression levels ..."
    #pappi.sql.dump_csv(hpa_levels_file, "hpa_gene_levels", con)

