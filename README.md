PAPPI - Protein Atlas and protein-protein interaction networks
==============================================================

The code consists of a python pipeline based on SQL (right now on SQLite in particular).

The pipeline imports all the various different data files, mostly in `comma separated values (CSV)` or
`tab separated values (TSV)` formats into a common SQL database. Then it performs
mapping of gene IDs from Entrez, Uniprot and Ensembl to HGNC symbols (names).
The HGNC symbols are used as main identifier in the further pipeline and analysis.

A few genes are still lost in this pipeline approach, as they can not be mapped to
any HGNC symbol.

The python pipeline is in the `src/` folder. The `init_pappi.py` script has to be called
in order to import all files. It is also there, that all the paths to the various
data files is configured.

Right now this script is still a bit messy and dependent from all the file paths (i.e. hardcoded filenames).


The `analysis` folder contains R scripts for analysis of the data.

Edit the sql_config.R file to point to the correct sqlite file.


Get the data:
----------

Create a folder for all the data files (preferably on an internal, and not a network drive).


### Human protein atlas

Go to the downloadable data section at [proteinatlas.org](http://www.proteinatlas.org/about/download)
and download the [normal_tissue.csv.zip](http://www.proteinatlas.org/download/normal_tissue.csv.zip) file.

Unpack it into the data folder.


### string-db

Go to [string-db](http://string-db.org/) to the Download section
and download the protein.links.vX.XX.txt.gz file. This file contains all protein protein
interactions using Ensembl IDs and a reliability score per interaction.

Direct download link (version 9.05): [here](http://string-db.org/newstring_download/protein.links.v9.05.txt.gz)


### CCSB network:

Sign up at the CCSB and download the Human interactome database [here](http://interactome.dfci.harvard.edu/H_sapiens/).

Save it into the data folder.


### MMC network

This is the protein complex network from the paper ["A Census of Human Soluble Protein Complexes"](http://www.cell.com/abstract/S0092-8674(12)01006-9)
by Havugimana et al. 

Download the network from the supplemental information Table S2 [here](http://www.sciencedirect.com/science/article/pii/S0092867412010069)
The protein-protein interactions are in the 

This network is imported as a CSV file holding only the two gene names, thus these have to be copied
into a separate file for import into the PAPPI pipeline.

Right now this network is already part of the repository. This might complicate giving free access to the
repo, as I am currently unsure about the licencing of the PPI network from Havugimana et al.



Mapping files
-------------

### BioMart Gene<->Protein mapping

Get the ensembl ID mapping from
http://www.ensembl.org/biomart/martview/
using homo sapiens version 71 database and output the attributes:

- Ensembl Gene ID
- Ensembl Protein ID

with NO filters


### BioMart Ensembl<->HGNC mapping
Go to: http://www.ensembl.org/biomart/martview

Choose "Ensembl Genes 71" and table "Homo sapiens genes"

Include following fields for the table:

- Ensembl Gene ID
- Associated Gene Name
- UniProt/SwissProt ID
- HGNC ID(s)
- EntrezGene ID

Export the table as CSV (and choose "Unique results only")


### HGNC mapping

A few files for Gene-ID mapping/matching need to be downloaded and imported as well.

Get data from: [http://www.genenames.org/cgi-bin/hgnc_stats]
Goto `Locus Group`: "protein-coding gene" and click "Custom".
Choose only the Columns:

- HGNC ID
- Approved Symbol
- Approved Name
- Status
- Entrez Gene ID
- Ensembl Gene ID

(and from external sources):

- Entrez Gene ID (supplied by NCBI)
- UniProt ID (supplied by UniProt)
- Ensembl ID (supplied by Ensembl)

Make sure to deselect (exclude) the status: "Entry and Symbol Withdrawn"

Full URL to results:
[BioMart Mapping](http://www.genenames.org/cgi-bin/hgnc_downloads?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=md_eg_id&col=md_prot_id&col=md_ensembl_id&status=Approved&status_opt=2&where=%28%28gd_pub_chrom_map+not+like+%27%25patch%25%27+and+gd_pub_chrom_map+not+like+%27%25ALT_REF%25%27%29+or+gd_pub_chrom_map+IS+NULL%29+and+gd_locus_type+%3D+%27gene+with+protein+product%27&order_by=gd_hgnc_id&format=text&limit=&submit=submit)

  

TODO:
-----

- write this readme, should include:
	+ download locations for all downloadable data (e.g. HPA, string-db, CCSB)
	+ how to import everything (which configs to edit how)
	+ overview of code layout (modules)
	+ explanation about different .R analysis scripts
	
- write a more detailed documentation !?
- ...

