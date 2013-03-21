
/* clear non-matching entires from the matching table */
DELETE FROM entrez_to_ensembl
WHERE EntrezID = '' OR EnsemblID = '';


/* 
 * create a ENSG based PPI network from the CCSB data
 * using the entrez_to_ensembl matching
 */
DROP TABLE IF EXISTS ppi_genes;
CREATE TABLE ppi_genes AS
SELECT
	b.EnsemblID AS Gene1,
	c.EnsemblID AS Gene2
FROM ccsb_raw AS a
INNER JOIN entrez_to_ensembl AS b
	ON a.Gene_IDA = b.EntrezID
INNER JOIN entrez_to_ensembl AS c
	ON a.Gene_IDB = c.EntrezID

