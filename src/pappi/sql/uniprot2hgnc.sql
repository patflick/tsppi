DROP TABLE IF EXISTS uniprot_to_hgnc_initial;
CREATE TABLE uniprot_to_hgnc_initial AS
SELECT 
	[HGNC.ID],
	[HGNC.Symbol],
	UniprotID
FROM hgnc_raw
WHERE UniprotID != '' AND UniprotID != '-'
;



/* 
 * now a problem remains, where one uniprot ID might map to multiple HGNC IDs/Symbols,
 * this is "solved" by mapping each uniprot ID to its minimum HGNC ID
 */

DROP TABLE IF EXISTS uniprot_to_hgnc_single;
CREATE TABLE uniprot_to_hgnc_single AS
SELECT
	MIN([HGNC.ID]) as [HGNC.ID],
	UniprotID
FROM uniprot_to_hgnc_initial
GROUP BY UniprotID
;

/* inner join, to select only those that are needed */
DROP TABLE IF EXISTS uniprot_to_hgnc;
CREATE TABLE uniprot_to_hgnc AS
SELECT
	a.[HGNC.ID] as [HGNC.ID],
	a.[HGNC.Symbol] as [HGNC.Symbol],
	b.UniprotID as UniprotID
FROM uniprot_to_hgnc_initial as a
INNER JOIN uniprot_to_hgnc_single as b
	ON a.[HGNC.ID] = b.[HGNC.ID] AND a.UniprotID = b.UniprotID
;

/* clean up, i.e. delete now unnecessary tables */
DROP TABLE uniprot_to_hgnc_initial;
DROP TABLE uniprot_to_hgnc_single;




