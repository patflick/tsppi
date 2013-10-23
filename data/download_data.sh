#!/bin/bash

# This script downloads all publicly available data into the correct folders.
# Run this before using the pipeline.
# Don't forget to download the data, that is not publicly available (see README.md for 
# sources and further explainations)

# A few settings
TMP_FOLDER=tmp


# make tmp folder
rm -R $TMP_FOLDER
mkdir -p $TMP_FOLDER


####################
#   PPI networks 
###################

# 1.) Bossi & Lehner:
BOSSI_URL=http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2683721/bin/msb200917-s2.zip
BOSSI_ZIP=$TMP_FOLDER/bossi-lehner-suppl.zip

curl -o $BOSSI_ZIP $BOSSI_URL
unzip $BOSSI_ZIP -d $TMP_FOLDER
cp $TMP_FOLDER/CRG-human-interactome/CRG.integrated.human.interactome.txt ppis/


