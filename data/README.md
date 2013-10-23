Data Sources
============

Protein-Protein Interaction Networks
------------------------------------

### Bossi & Lehner composite PPI network

From the supplementary section of the paper ["Tissue specificity and the human protein interaction network"](http://www.ncbi.nlm.nih.gov/pubmed/19357639).

This can be automatically downloaded via the `download_data.sh` script.



### Dana Farber CCSB PPI network

The Human Interactome 2012 is still in preliminary form, thus you have
to sign up at the CCSB and to download the Human interactome database [here](http://interactome.dfci.harvard.edu/H_sapiens/).

Download the `HI_2012_PRE.tsv` file and save it into the `ppis/` folder.



### Havugimana et al. protein complexes PPI network

This is the protein complex network from the paper ["A Census of Human Soluble Protein Complexes"](http://www.ncbi.nlm.nih.gov/pubmed/22939629)
by Havugimana et al. 

Download the network from the supplemental information Table S2 [here](http://www.sciencedirect.com/science/article/pii/S0092867412010069)

The protein-protein interactions are in the Excel sheet "14K Denoised PPI".
Copy the first two columns into a tab-separated text file with headers
"Gene1" and "Gene2".

#### Note:

Right now this network is already part of the repository. This might complicate giving free access to the
repo, as I am currently unsure about the licencing of the PPI network from Havugimana et al.



### string-db

The protein interaction data can be loaded from string-db
automatically with the `download_data.sh` shell script.

If you already downloaded the string-db file or import
it manually for other reasons, do the following:

Go to [string-db](http://string-db.org/) to the Download section
and download the protein.links.vX.XX.txt.gz file. This file contains all protein protein
interactions using Ensembl IDs and a reliability score per interaction.

Direct download link (for version 9.05): [here](http://string-db.org/newstring_download/protein.links.v9.05.txt.gz)

Save this file into the `data/download` folder, and then run
the `download_data.sh` script, which will pre-filter the
string-db interactions for human interactions with score >= 0.7




Protein Expression data sets
----------------------------

These are all downloaded and unpacked automatically via the
`download_data.sh` script. The following just describes the
data sources.


### Human Protein Atlas (HPA)

From the [Human Protein Atlas](http://www.proteinatlas.org/) version 11 
the `normal_tissue.csv` and `subcellular_location.csv` files are used.

Manual download [here](http://www.proteinatlas.org/about/download)


### GeneAtlas (now BioGPS)

This is RNA micro-array data from Su Al et al. 2009.
It can be downloaded from BioGPS [here](http://biogps.org/downloads/).

The files `gnf1h-gcrma.zip` and `gnf1h-anntable.zip` are needed, where
the first is the actual expression data while the second holds the annotation
for the genes.


### Illumina Body Map RNAseq

The RNAseq data from Illumina Body Map can be downloaded from EBI
[here](http://www-test.ebi.ac.uk/gxa/experiments/E-MTAB-513).

The automatic script removes all filters (no gene selection, especially not only protein coding genes; and a cutoff of 0).


