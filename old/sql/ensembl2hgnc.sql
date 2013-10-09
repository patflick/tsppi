DROP TABLE IF EXISTS ensembl_to_hgnc_1;
CREATE TABLE ensembl_to_hgnc_1 AS
SELECT 
	[HGNC.ID],
	[HGNC.Symbol],
	CASE
		WHEN [Ensembl.EnsemblID] = '' THEN [HGNC.EnsemblID]
		ELSE [Ensembl.EnsemblID]
	END AS EnsemblID
FROM hgnc_raw
/* this looses 19060 - 18895 genes */
WHERE EnsemblID != ''
;


/* 165 genes in HGNC don't have a match in Ensembl */
/* FIXME this doesn't belong here, this comes into hgnc_to_ensembl, not the other way around (as here) */
/*
DROP TABLE IF EXISTS hgnc_no_ensembl;
CREATE TABLE hgnc_no_ensembl AS
SELECT
	[HGNC.ID],
	[HGNC.Symbol]
FROM
	hgnc_raw
	WHERE [Ensembl.EnsemblID] = '' AND [HGNC.EnsemblID] = ''
;
*/

/* get biomart entries that are not in HGNC */
DROP TABLE IF EXISTS ensembl_not_in_hgnc;
CREATE TABLE ensembl_not_in_hgnc AS
SELECT DISTINCT
	Ensembl_GeneID,
	Gene_Name
FROM biomart_raw
WHERE Ensembl_GeneID NOT IN
(
	SELECT EnsemblID FROM ensembl_to_hgnc_1
)
;

/* now try to join on Gene Name to hgnc, to get synonym gene names (ensembl synonyms) into the mapping table */
DROP TABLE IF EXISTS ensembl_to_hgnc_2;
CREATE TABLE ensembl_to_hgnc_2 AS
SELECT
	a.[HGNC.ID] as [HGNC.ID],
	a.[HGNC.Symbol] as [HGNC.Symbol],
	b.Ensembl_GeneID as EnsemblID
FROM hgnc_raw as a
INNER JOIN
	ensembl_not_in_hgnc as b
ON
	a.[HGNC.Symbol] = b.Gene_Name
;


/* union tables to get an ensembl to hgnc mapping */
DROP TABLE IF EXISTS ensembl_to_hgnc_3;
CREATE TABLE ensembl_to_hgnc_3 AS
SELECT * FROM ensembl_to_hgnc_1
UNION
SELECT * FROM ensembl_to_hgnc_2
;



/*  REMOVE DUPLICATES! */

/* I now have multiple matches per Ensembl ID in HGNC 
 * mostly pairs: Name-1 and Name-2 or Name-A and Name-B
 */
DROP TABLE IF EXISTS duplicates_ensembl;
CREATE TABLE duplicates_ensembl AS
select * FROM ensembl_to_hgnc_3
WHERE EnsemblID IN
(
select EnsemblID FROM ensembl_to_hgnc_3
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


/* Create final table! */

DROP TABLE IF EXISTS ensembl_to_hgnc;
CREATE TABLE ensembl_to_hgnc AS
SELECT * FROM ensembl_to_hgnc_3
EXCEPT
SELECT * FROM duplicates_delete


