"""
Defines the file paths for all data files
"""

# for getting the current file path
import os

# go two folders up
cur_folder = os.path.dirname(__file__)
base_folder = os.path.join(cur_folder, '..', '..')

EXT_DATA_FOLDER = os.path.join(base_folder, '..', 'data')
DATA_FOLDER = os.path.join(base_folder, 'data')

PPI_DATA_FOLDER = os.path.join(DATA_FOLDER, 'ppis')
EXPR_DATA_FOLDER = os.path.join(DATA_FOLDER, 'expr')

# TODO add these files into the repo data folder !?
HGNC_FILE = os.path.join(EXT_DATA_FOLDER, 'hgnc_downloads.txt')
BIOMART_FILE = os.path.join(EXT_DATA_FOLDER, 'mart_export.csv')
DATABASE = os.path.join(EXT_DATA_FOLDER, 'test_matching.sqlite')

GNF1H_ANNOT_FILE = os.path.join(EXPR_DATA_FOLDER, 'gnf1h.annot2007.tsv')
U133A_ANNOT_FILE = os.path.join(EXPR_DATA_FOLDER, 'GPL96-15653.txt')

# paths to PPIs
CCSB_FILE = os.path.join(PPI_DATA_FOLDER, 'HI_2012_PRE.tsv')
BOSSI_FILE = os.path.join(PPI_DATA_FOLDER, 'CRG.integrated.human.interactome.txt')
HAVU_FILE = os.path.join(PPI_DATA_FOLDER, 'cell_havugimana_ppi.tsv')

# import expression data sets
HPA_FILE = os.path.join(EXPR_DATA_FOLDER, 'normal_tissue.csv')
EMTAB_FILE = os.path.join(EXPR_DATA_FOLDER, 'E-MTAB-513.tsv')
RNASEQ_ATLAS_FILE = os.path.join(EXPR_DATA_FOLDER, 'RNA_Seq_Atlas_rev1.txt')
GENE_ATLAS_FILE = os.path.join(EXPR_DATA_FOLDER, 'U133AGNF1B.gcrma.avg.csv')
