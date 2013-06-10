
DROP TABLE IF EXISTS entrez_to_hgnc;
CREATE TABLE entrez_to_hgnc AS
SELECT
	[HGNC.ID],
	[HGNC.Symbol],
    CASE
		WHEN [NCBI.EntrezID] = '' THEN [HGNC.EntrezID]
		ELSE [NCBI.EntrezID]
	END AS EntrezID
FROM hgnc_raw
WHERE EntrezID != ''

/* TODO throw out duplicates */




