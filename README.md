PAPPI - Protein Atlas and protein-protein interaction networks
==============================================================

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





  

TODO:
-----

- write this readme, should include:
	+ download locations for all downloadable data (e.g. HPA, string-db, CCSB)
	+ how to import everything (which configs to edit how)
	+ overview of code layout (modules)
	+ explanation about different .R analysis scripts
	
- write a more detailed documentation !?
- ...

