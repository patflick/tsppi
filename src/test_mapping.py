import pappi.id_mapping
import pappi.sql
import os

DATA_FOLDER = '/home/patrick/dev/bio/data/'

HGNC_FILE = DATA_FOLDER + 'hgnc_entrez_ensembl_uniprot.txt'
BIOMART_FILE = DATA_FOLDER + 'mart_export.csv'
DATABASE = DATA_FOLDER + 'test_matching.sqlite'


# first delete an old DB, to make sure everything is new
if (os.path.exists(DATABASE)):
    os.remove(DATABASE)

# get new database connection
con = pappi.sql.get_conn(DATABASE)
pappi.id_mapping.import_biomart_file(BIOMART_FILE, con)
pappi.id_mapping.import_hgnc_file(HGNC_FILE, con)
pappi.id_mapping.create_mapping_table("ensembl", "hgnc", con, True)
pappi.id_mapping.create_mapping_table("uniprot", "hgnc", con, True)
pappi.id_mapping.create_mapping_table("entrez", "hgnc", con, True)
pappi.id_mapping.create_mapping_table("entrez", "uniprot", con, True)
