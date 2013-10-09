

import pappi

DATABASE='/cygdrive/d/PPI/hpaDB.sqlite'

PPI_FILE='/cygdrive/d/PPI/string-db/human_protein.links.v9.0_700.csv'
HPA_FILE='../../hpa/data/normal_tissue.csv'
P2G_FILE='../../hpa/data/ensembl_ID_matching.csv'

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
pappi.import_hpa(HPA_FILE, DATABASE)
pappi.import_ppi(PPI_FILE, P2G_FILE, DATABASE)