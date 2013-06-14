CREATE TABLE hgnc_no_ensembl AS
SELECT * FROM hgnc_raw WHERE [HGNC.EnsemblID] = "" AND [Ensembl.EnsemblID] = "";


CREATE TABLE test_blah AS
SELECT distinct a.[HGNC.ID] as HGNC_ID, b.[Ensembl_GeneID] as EnsemblID FROM hgnc_no_ensembl AS a INNER JOIN biomart_raw AS b ON a.[HGNC.ID] = b.HGNC_ID;

SELECT distinct a.[HGNC.ID] as HGNC_ID, a.[HGNC.Symbol] as Name, b.[Ensembl_GeneID] as EnsemblID FROM hgnc_no_ensembl AS a INNER JOIN biomart_raw AS b ON a.[HGNC.Symbol] = b.Gene_Name



SELECT * FROM hgnc_no_ensembl WHERE [NCBI.EntrezID] IN (SELECT DISTINCT Gene_IDA FROM ccsb_raw)

/* HPA genes NOT in the HGNC database, returns 565 !!! */
SELECT distinct Gene FROM hpa_normal_tissue WHERE Gene NOT IN (SELECT [Ensembl.EnsemblID] FROM hgnc_raw UNION SELECT [HGNC.EnsemblID] FROM hgnc_raw)


/* HPA genes NOT in the biomart results
 * with all filters:     396 !!! 
 * only transcript >= 1: 47
 *   but those are mostly: uncharacterized proteins (no external references to entrez, hgnc or uniprot)
 * 
 */
SELECT distinct Gene FROM hpa_normal_tissue WHERE Gene NOT IN (SELECT Ensembl_GeneID FROM biomart_raw)




/*
 * Entrez Genes in CCSB PPI that are not in the biomart table
 * 47 Genes, leading to 448 interactions
 
 */

SELECT distinct Gene_IDA FROM ccsb_raw WHERE Gene_IDA NOT IN (SELECT EntrezID FROM biomart_raw)
SELECT Gene_IDA, Gene_IDB FROM ccsb_raw WHERE Gene_IDA NOT IN (SELECT EntrezID FROM biomart_raw) OR Gene_IDB NOT IN (SELECT EntrezID FROM biomart_raw) 


/* 
 * same with hgnc:
 * total of 116 Genes, involved in 412 interactions!
 * out of the 116 genes:
 *    47 pseudogenes
 *    18 "non-protein coding"
 *    10 antisense RNA
 *    10 uncharacterized
 *    
 */
SELECT DISTINCT Gene
FROM
(
SELECT Gene_IDA as Gene FROM ccsb_raw WHERE Gene_IDA NOT IN (SELECT [NCBI.EntrezID] FROM hgnc_raw)
UNION
SELECT Gene_IDB as Gene FROM ccsb_raw WHERE Gene_IDB NOT IN (SELECT [NCBI.EntrezID] FROM hgnc_raw)
)


/* CCSB genes in HGNC that don't have an Ensembl ID:
 *  ->  24
 */
SELECT * FROM hgnc_raw WHERE [NCBI.EntrezID] IN
(
SELECT Gene_IDA as Gene FROM ccsb_raw
UNION
SELECT Gene_IDB as Gene FROM ccsb_raw
)
AND [Ensembl.EnsemblID] = ""
AND [HGNC.EnsemblID] = ""




/* check the biomart for the 565 genes from HPA that are not in HGNC ensembl */

SELECT * FROM biomart_raw WHERE Ensembl_GeneID IN
(
SELECT distinct Gene FROM hpa_normal_tissue WHERE Gene NOT IN (SELECT [Ensembl.EnsemblID] FROM hgnc_raw UNION SELECT [HGNC.EnsemblID] FROM hgnc_raw)
)


SELECT * FROM hgnc_raw WHERE [HGNC.Symbol] = "PNRC2"



/* check/verify my new table */

/* still 263 Ensembl genes lost */
SELECT DISTINCT Gene, Reliability FROM hpa_normal_tissue WHERE Gene NOT IN (SELECT EnsemblID FROM ensembl_to_hgnc)

/* but 226 are in biomart! */
SELECT DISTINCT Ensembl_GeneID FROM biomart_raw WHERE Ensembl_GeneID IN (
SELECT DISTINCT Gene FROM hpa_normal_tissue WHERE Gene NOT IN (SELECT EnsemblID FROM ensembl_to_hgnc)
)

/* those: */
SELECT * FROM biomart_raw WHERE Ensembl_GeneID IN (
SELECT DISTINCT Gene FROM hpa_normal_tissue WHERE Gene NOT IN (SELECT EnsemblID FROM ensembl_to_hgnc)
)

/* get entrez ids of those
 *  this can only get me another 8 HGNC entries :(
 */
CREATE TABLE tmp_entrez_biomart AS
SELECT DISTINCT EntrezID FROM biomart_raw WHERE Ensembl_GeneID IN (
SELECT DISTINCT Gene FROM hpa_normal_tissue WHERE Gene NOT IN (SELECT EnsemblID FROM ensembl_to_hgnc)
)
AND EntrezID != ""




/* I now have multiple matches per Ensembl ID in HGNC 
 * mostly pairs: Name-1 and Name-2 or Name-A and Name-B
 */
DROP TABLE IF EXISTS duplicates_ensembl;
CREATE TABLE duplicates_ensembl AS
select * FROM ensembl_to_hgnc
WHERE EnsemblID IN
(
select EnsemblID FROM ensembl_to_hgnc
GROUP BY EnsemblID
HAVING count() >= 2
)
ORDER BY EnsemblID
;

/* for each duplicate, keep the one with the minimum hgnc id */
DROP TABLE IF EXISTS duplicates_keep_tmp;
CREATE TABLE duplicates_keep_tmp AS
SELECT MIN([HGNC.ID]) as [HGNC.ID], EnsemblID FROM duplicates_ensembl
GROUP BY EnsemblID
;

DROP TABLE IF EXISTS duplicates_keep;
CREATE TABLE duplicates_keep AS
SELECT a.* FROM duplicates_ensembl AS a
INNER JOIN duplicates_keep_tmp AS b
ON a.[HGNC.ID] = b.[HGNC.ID] AND a.EnsemblID = b.EnsemblID
;

/* create table of duplicates that have to be deleted from the main table */
DROP TABLE IF EXISTS duplicates_delete;
CREATE TABLE duplicates_delete AS
SELECT * FROM duplicates_ensembl
EXCEPT
SELECT * FROM duplicates_keep
;




select * FROM ensembl_to_hgnc
GROUP BY EnsemblID
HAVING count() >= 2
ORDER BY EnsemblID



/* check for multiples in uniprot2hgnc mapping */
select * FROM uniprot_to_hgnc
WHERE UniprotID IN
(
select UniprotID FROM uniprot_to_hgnc
GROUP BY UniprotID 
HAVING count() >= 2
)
ORDER BY UniprotID 




SELECT * FROM hpa_gene_levels WHERE CountTotal > 80

