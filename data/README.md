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