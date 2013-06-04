'''
Configuration for the PAPPI package.
Sets the paths and filenames for the SQL scripts in the ./sql folder.

@author: Patrick Flick
'''

import os

# the SQL folder path, relative to this file
PAPPI_SQL_FOLDER_PATH=os.path.abspath(os.path.join(os.path.dirname(__file__), 'sql'))

# HPA import, initialization and processing
PAPPI_SQL_HPA_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'hpa_filter.sql')
PAPPI_SQL_HGNC_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'hgnc_filter.sql')
PAPPI_SQL_HPA_GENE_LEVELS_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'hpa_gene_levels.sql')

# PPI import and initialization
PAPPI_SQL_STRINGDB_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'stringdb_filter.sql')
PAPPI_SQL_CCSB_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'ccsb_filter.sql')
PAPPI_SQL_MMC_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'ppi_mmc_filter.sql')

# PPI processing
PAPPI_SQL_EDGE_EXPR_SCRIPT =  os.path.join(PAPPI_SQL_FOLDER_PATH,'edge_expression.sql')

# mapping initialization scripts (From other IDs to HGNC ID)
PAPPI_SQL_ENSEMBL2HGNC_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'ensembl2hgnc.sql')
PAPPI_SQL_UNIPROT2HGNC_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'uniprot2hgnc.sql')

# default SQLite database file
PAPPI_SQLITE_DEFAULT_DB = 'pappiDB.sqlite'