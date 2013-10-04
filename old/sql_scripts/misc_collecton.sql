
/* create table of all ENSG reliable genes from the HPA */
DROP TABLE IF EXISTS hpa_reliable_ensg_genes;
CREATE TABLE hpa_reliable_ensg_genes AS

SELECT Gene FROM hpa_normal_tissue_raw WHERE 
    
	(
	Reliability = "High"
        OR
        Reliability = "Medium"
        OR
        Reliability = "Supportive"
    )
GROUP BY Gene;

/* 
 * merge the previous table with the HGNC Symbols,
 * resulting in a list of all reliable genes as (ENSG, HGNC_Symbol).
 */
DROP TABLE IF EXISTS hpa_reliable_hgnc_and_ensg;
CREATE TABLE hpa_reliable_hgnc_and_ensg AS
SELECT
    a.Gene as ENSG_Gene,
    b.[HGNC.Symbol] as HGNC_Symbol
FROM hpa_reliable_ensg_genes as a
INNER JOIN ensembl_to_hgnc as b
ON a.Gene = b.EnsemblID
;
