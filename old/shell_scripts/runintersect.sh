#!/bin/bash

INTERSECT_PY=../python/src/intersect_hpa.py
PPI_FILTER_PY=../python/src/ppi_filter_genes.py
PPI_GENES=/cygdrive/d/PPI/string-db/human_protein.links.v9.0_700_genes.csv

FILE1=../hpa/data/tissues/cerebral.cortex_neuronal.cells.csv
FILE2=../hpa/data/tissues/liver_hepatocytes.csv



OUTFOLDER=${FILE1%/*}/intersected

OUTFILE1=$OUTFOLDER/${FILE1##*/}
OUTFILE2=$OUTFOLDER/${FILE2##*/}
OUT_PPI=$OUTFOLDER/ppi.csv

mkdir $OUTFOLDER

# intersecting tissue datasets
echo "intersecting tissue datasets..."
python $INTERSECT_PY -o $OUTFOLDER/genes.csv $FILE1 $FILE2 $OUTFILE1 $OUTFILE2

# intersecting resulting genes with PPI network
echo "intersecting resulting genes with PPI network..."
python $PPI_FILTER_PY $PPI_GENES $OUTFOLDER/genes.csv $OUT_PPI

# done
echo "DONE!"