#!/bin/bash

# Script Variables: EDIT THOSE:

# minium score for edges in the PPI
MIN_SCORE=700

# the input PPI file:
STRING_DB_IN=/cygdrive/d/PPI/string-db/protein.links.v9.0.txt.gz

# the input ENSG to ENSP mapping file:
ENSG_ENSP_MAPPING=$HOME/dev/ppi/hpa/data/ensemblID_matching_cleared.csv

# technical settings:

# the python script to collapse the PPI from proteins -> genes
P2G_MAPPING_PY=../python/src/ppi_P2G.py

# DO NOT EDIT BELOW THIS LINE:

# extract folder and file
PPI_FOLDER=${STRING_DB_IN%/*}
PPI_FILE=${STRING_DB_IN##*/}
# set output file for the human PPI
HUMAN_PPI=$PPI_FOLDER/human_${PPI_FILE%.txt.gz}_$MIN_SCORE.csv

# filter PPI for only human and a given minimum score
if [ ! -f $HUMAN_PPI ];
then
	echo "Filtering PPI for HUMAN and SCORE >= $MIN_SCORE ..."
	zgrep ^"9606\." $STRING_DB_IN | awk -v score=$MIN_SCORE '($3 >= score) {print}' > $HUMAN_PPI
else
	echo "File ${HUMAN_PPI##*/} already exists, skipping the first filtering step."
fi

# replace proteins->genes and filter duplicates using the python script ppi_P2G:
HUMAN_PPI_GENES=${HUMAN_PPI%.csv}_genes.csv
echo "replacing ENSP -> ENSG and removing duplicates ..."
python $P2G_MAPPING_PY $HUMAN_PPI $ENSG_ENSP_MAPPING $HUMAN_PPI_GENES