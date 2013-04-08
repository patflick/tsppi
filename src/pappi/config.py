'''
Configuration for the PAPPI package.

@author: Patrick Flick
'''

import os

# the SQL folder path, relative to this file
PAPPI_SQL_FOLDER_PATH=os.path.abspath(os.path.join(os.path.dirname(__file__), 'sql'))
PAPPI_SQL_STRINGDB_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'stringdb_filter.sql')
PAPPI_SQL_CCSB_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'ccsb_filter.sql')
PAPPI_SQL_HPA_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'hpa_filter.sql')
PAPPI_SQL_HGNC_FILTER_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'hgnc_filter.sql')
PAPPI_SQL_HPA_GENE_LEVELS_SCRIPT = os.path.join(PAPPI_SQL_FOLDER_PATH,'hpa_gene_levels.sql')
PAPPI_SQLITE_DEFAULT_DB = 'pappiDB.sqlite'