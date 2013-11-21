#!/bin/bash

# This script downloads all publicly available data into the correct folders.
# Run this before using the pipeline.
# Don't forget to download the data, that is not publicly available (see README.md for 
# sources and further explainations)

# A few settings
TMP_FOLDER=tmp
DOWNLOAD_FOLDER=download

# make tmp folder
rm -R $TMP_FOLDER
mkdir -p $TMP_FOLDER
mkdir -p $DOWNLOAD_FOLDER

mkdir ppis
mkdir expr

####################
#   PPI networks 
###################

# 1.) Bossi & Lehner:
BOSSI_URL=http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2683721/bin/msb200917-s2.zip
BOSSI_ZIP=$DOWNLOAD_FOLDER/bossi-lehner-suppl.zip

if [ -f $BOSSI_ZIP ]
then
    echo "Bossi & Lehner file already exists, skipping download"
else
    echo "Downloading Bossi & Lehner PPI network"
    curl -o $BOSSI_ZIP $BOSSI_URL
fi

unzip $BOSSI_ZIP -d $TMP_FOLDER
cp $TMP_FOLDER/CRG-human-interactome/CRG.integrated.human.interactome.txt ppis/

# 2.) CCSB (only HI-2011, the HI-2012-Pre is not yet available)

# TODO when HI-2012 becomes publicly available
#      at the moment the HI-2012-Pre has to be downloaded manually (requires
#      an account to download)

# 3.) string-db
STRING_VERSION=9.05
STRING_CUTOFF=500 # cutoff at a reliability of 700 for import
STRING_URL=http://string-db.com/newstring_download/protein.links.v$STRING_VERSION.txt.gz
STRING_GZ=$DOWNLOAD_FOLDER/protein.links.v$STRING_VERSION.txt.gz
STRING_FILE=./ppis/protein.links.v$STRING_VERSION.$STRING_CUTOFF.txt

if [ -f $STRING_GZ ]
then
    echo "string-db file already exists, skipping download"
else
    echo "Downloading string-db file"
    curl -o $STRING_GZ $STRING_URL
fi

# pre-filter string-db from the gz into an importable text file
if [ ! -f $STRING_FILE ]
then
    zgrep ^"9606\." $STRING_GZ | awk -v score=$STRING_CUTOFF '($3 >= score) {print}' > $STRING_FILE
fi


########################
# expression data sets
########################

# 1.) Human Protein Atlas data
# expression data
HPA_URL=http://www.proteinatlas.org/download/normal_tissue.csv.zip
HPA_ZIP=normal_tissue.csv.zip
HPA_FILE=normal_tissue.csv

# localization data
HPA_LOC_URL=http://www.proteinatlas.org/download/subcellular_location.csv.zip
HPA_LOC_ZIP=subcellular_location.csv.zip
HPA_LOC_FILE=subcellular_location.csv

if [ -f $DOWNLOAD_FOLDER/$HPA_ZIP ]
then
    echo "HPA normal_tissue file already exists, skipping download"
else
    echo "Downloading HPA normal_tissue data"
    curl -o $DOWNLOAD_FOLDER/$HPA_ZIP $HPA_URL
fi


if [ -f $DOWNLOAD_FOLDER/$HPA_LOC_ZIP ]
then
    echo "HPA subcellular_location file already exists, skipping download"
else
    echo "Downloading HPA subcellular_location data"
    curl -o $DOWNLOAD_FOLDER/$HPA_LOC_ZIP $HPA_LOC_URL
fi

# unzip and move the files to the expr folder
if [ ! -f expr/$HPA_FILE ]
then
    unzip $DOWNLOAD_FOLDER/$HPA_ZIP -d $TMP_FOLDER
    mv $TMP_FOLDER/$HPA_FILE expr/
fi

if [ ! -f expr/$HPA_LOC_FILE ]
then
    unzip $DOWNLOAD_FOLDER/$HPA_LOC_ZIP -d $TMP_FOLDER
    mv $TMP_FOLDER/$HPA_LOC_FILE expr/
fi


# 2.) Gene Atlas (Su Al et al. 2009) from BioGPS
# expression file
GNF1H_URL=http://plugins.biogps.org/download/gnf1h-gcrma.zip
GNF1H_ZIP=gnf1h-gcrma.zip
GNF1H_FILE=U133AGNF1B.gcrma.avg.csv
# gene annotations
GNF1H_ANN_URL=http://plugins.biogps.org/download/gnf1h-anntable.zip
GNF1H_ANN_ZIP=gnf1h-anntable.zip
GNF1H_ANN_FILE=gnf1h.annot2007.tsv
# more annotations from U133A
U133A_ANN_URL="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL96&id=15653&db=GeoDb_blob82"
U133A_ANN_FILE=GPL96-15653.txt

if [ -f $DOWNLOAD_FOLDER/$GNF1H_ZIP ]
then
    echo "GeneAtlas GNF1H file already exists, skipping download"
else
    echo "Downloading GeneAtlas GNF1H file"
    curl -o $DOWNLOAD_FOLDER/$GNF1H_ZIP $GNF1H_URL
fi

if [ -f $DOWNLOAD_FOLDER/$GNF1H_ANN_ZIP ]
then
    echo "GeneAtlas GNF1H annotations file already exists, skipping download"
else
    echo "Downloading GeneAtlas GNF1H annotations file"
    curl -o $DOWNLOAD_FOLDER/$GNF1H_ANN_ZIP $GNF1H_ANN_URL
fi

# unpack files and move/copy to data/expr folder
if [ ! -f expr/$GNF1H_FILE ]
then
    unzip $DOWNLOAD_FOLDER/$GNF1H_ZIP -d $TMP_FOLDER
    # prepend the file with the column name for the first column (which is empty)
    echo -n "Gene_ID" | cat - $TMP_FOLDER/$GNF1H_FILE > ./expr/$GNF1H_FILE
fi

if [ ! -f expr/$GNF1H_ANN_FILE ]
then
    unzip $DOWNLOAD_FOLDER/$GNF1H_ANN_ZIP -d $TMP_FOLDER
    mv $TMP_FOLDER/$GNF1H_ANN_FILE expr/
fi

# download the U133A annotation file
if [ ! -f ./expr/$U133A_ANN_FILE ]
then
    curl -o expr/$U133A_ANN_FILE $U133A_ANN_URL
fi


# 3.) RNAseq Illumina Body Map from EBI

EBI_MTAB_URL="http://www-test.ebi.ac.uk/gxa/experiments/E-MTAB-513.tsv?cutoff=0&geneQuery="
EBI_MTAB_FILE=expr/E-MTAB-513.tsv

if [ -f $EBI_MTAB_FILE ]
then
    echo "EBI MTAB file already exists, skipping download"
else
    echo "Downloading EBI E-MTAB-513 file"
    curl -o $EBI_MTAB_FILE $EBI_MTAB_URL
fi


# 4.) RNA Seq Atlas from medicalgenomics.org
RNASEQ_ATLAS_URL="http://medicalgenomics.org/rna_seq_atlas/download?download_revision1=1"
RNASEQ_ATLAS_FILE=RNA_Seq_Atlas_rev1.txt

if [ -f expr/$RNASEQ_ATLAS_FILE ]
then
    echo "RNA Seq Atlas file already exists, skipping download"
else
    echo "Downloading RNA Seq Atlas file"
    curl -o expr/$RNASEQ_ATLAS_FILE $RNASEQ_ATLAS_URL
fi

