

from pappi import hpa
from pappi import ppi
import os

DATABASE='/cygdrive/d/PPI/hpaDB.sqlite'

PPI_FILE='/cygdrive/d/PPI/string-db/human_protein.links.v9.0_700.csv'
HPA_FILE='/home/flick/dev/ppi/hpa/data/normal_tissue.csv'
P2G_FILE='/home/flick/dev/ppi/hpa/data/ensembl_ID_matching.csv'

# TODO properly test this
# TODO continue with pipeline (merge two tissues, intersect with PPI, do graph scoring)
# TODO actually i can do graph scoring for _every_ tissue as pre-work
# and then comparing is as simple as loading and diffing two tables in the DB 
# TODO read the HPA paper (look how they score the difference between tissues)
# TODO differences in interaction?? (
#        Tissue1 has protein A
#        Tissue2 as proteins A and B
#        Tissue3 has protein B
#        PPI has edge A<->B
#    => Proteins are not tissue specific, but interaction is tissue specific!
#    => How can I test for that?

# first delete an old DB, to make sure everything is new
os.remove(DATABASE)

hpa_file = open(HPA_FILE)
hpa.import_hpa(hpa_file, DATABASE)
ppi_file = open(PPI_FILE)
p2g_file = open(P2G_FILE)
ppi.import_stringdb(ppi_file, p2g_file, DATABASE)