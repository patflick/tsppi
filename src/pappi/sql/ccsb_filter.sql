
/* 
 * create a ENSG based PPI network from the CCSB data
 * using the entrez_to_ensembl matching
 */
DROP TABLE IF EXISTS ppi_genes;
CREATE TABLE ppi_ccsb_hgnc AS
SELECT
	b.[HGNC.Symbol] AS Gene1,
	c.[HGNC.Symbol] AS Gene2
FROM ccsb_raw AS a
INNER JOIN entrez_to_hgnc AS b
	ON a.Gene_IDA = b.EntrezID
INNER JOIN entrez_to_hgnc AS c
	ON a.Gene_IDB = c.EntrezID
;


DROP TABLE IF EXISTS ppi_hgnc;
CREATE TABLE ppi_hgnc AS
SELECT * FROM ppi_ccsb_hgnc
;
