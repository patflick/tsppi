

/* map the UniprotIDs of the MMC network to HGNC Symbols */

DROP TABLE IF EXISTS ppi_mmc_hgnc;
CREATE TABLE ppi_mmc_hgnc AS
SELECT
	b.[HGNC.Symbol] as Gene1,
	c.[HGNC.Symbol] as Gene2
FROM mmc_raw as a
INNER JOIN uniprot_to_hgnc as b
	ON a.Gene1 = b.UniprotID
INNER JOIN uniprot_to_hgnc as c
	ON a.Gene2 = c.UniprotID
;

/* activate mmc as default PPI */
DROP TABLE IF EXISTS ppi_hgnc;
CREATE TABLE ppi_hgnc AS
SELECT * FROM ppi_mmc_hgnc;

