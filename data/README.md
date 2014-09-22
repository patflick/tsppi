Data Sources
============

The majority of data sources are either provided in this repository or
automatically downloaded via scripts. However, some need to be manually obtained
and saved into the folder structure prior to executing the analysis pipeline.

### Automatic download

To download most of the data sources, simply execute:

```sh
./download.sh
```

and (**note:** requires python 3 to run):

```sh
python3 download_psicquic.py
```

### Manual download

The following data sources need to be manually obtained, since they are not
(yet) publicly available:

- Dana Farber CCSB HI-2012 PPI network (see below)


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

Right now this network is already part of the repository. This might complicate
giving free access to the repo, as I am currently unsure about the licensing of
the PPI network from Havugimana et al.


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

### PSICQUIC Composite network

Run the `download_psicquic.py` python script in this folder in order to download
the different *PSICQUIC* provided PPI networks.


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

### RNAseq Atlas from medicalgenomics.org

This data is also publicly available and can be downloaded from
[here](http://medicalgenomics.org/rna_seq_atlas).

This is also automatically downloaded with the `download.sh` script.


Mapping files
-------------

### BioMart id mapping
Go to: [ensembl.org](http://www.ensembl.org/biomart/martview)

Choose "Ensembl Genes 71" (or current version) and table "Homo sapiens genes"

Include following fields for the table:

- Ensembl Gene ID
- Ensembl Protein ID
- Associated Gene Name
- UniProt/SwissProt ID
- HGNC ID(s)
- EntrezGene ID

Export the table as CSV (and choose "Unique results only") and save this into
the file `mapping/mart_export.csv`.

Currently BioMart version 71 is provided in the repository.

### HGNC mapping

Get data from: [genenames.org](http://www.genenames.org/cgi-bin/hgnc_stats)
Goto `Locus Group`: "protein-coding gene" and click "Custom".
Choose only the Columns:

- HGNC ID
- Approved Symbol
- Approved Name
- Status
- Entrez Gene ID
- Ensembl Gene ID
(and from external sources)
- Entrez Gene ID (supplied by NCBI)
- UniProt ID (supplied by UniProt)
- Ensembl ID (supplied by Ensembl)

Make sure to deselect (exclude) the status: "Entry and Symbol Withdrawn"


Full URL to results:
[HGNC Mapping](http://www.genenames.org/cgi-bin/hgnc_downloads?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=md_eg_id&col=md_prot_id&col=md_ensembl_id&status=Approved&status_opt=2&where=%28%28gd_pub_chrom_map+not+like+%27%25patch%25%27+and+gd_pub_chrom_map+not+like+%27%25ALT_REF%25%27%29+or+gd_pub_chrom_map+IS+NULL%29+and+gd_locus_type+%3D+%27gene+with+protein+product%27&order_by=gd_hgnc_id&format=text&limit=&submit=submit)

This should return all 19060 genes.

Save this file into the `mapping/hgnc_downloads.txt`.
A version of this file is provided in the repository (not guaranteed to be up
to date).
