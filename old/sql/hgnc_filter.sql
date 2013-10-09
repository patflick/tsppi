

/*

Choosing the "right" Entrez ID
----------------------------

Running the following query on the raw data, gives rise to 3 mismatches.
In all 3 cases the NCBI data seems to be more up-to-date.
Thus in the following script, in case of mismatch, the NCBI ID is choosen
and the HGNC supplied EntrezID is discarded.

	SELECT * FROM hgnc_raw
	WHERE [HGNC.EntrezID] != [NCBI.EntrezID]
	AND [HGNC.EntrezID] != '' AND [NCBI.EntrezID] != ''

*/


/*

Choosing the "right" Ensembl ID
----------------------------

Running the following query on the raw data, returns 8 mismatches.
In all 8 cases the error lies with the HGNC data, thus in case
of a mismatch (same as for Entrez) the Ensembl data source is used
and the HGNC EnsemlbID discarded.

	SELECT * FROM hgnc_raw
	WHERE
		[HGNC.EnsemblID] != [Ensembl.EnsemblID]
	AND [HGNC.EnsemblID] != ''
	AND [Ensembl.EnsemblID] != ''

*/


/*
 * Create the Entrez to Ensembl matching. Always use the external (NCBI or Ensembl)
 * source, unless this is = '' (i.e. empty). 
 */

DROP TABLE IF EXISTS entrez_to_ensembl;
CREATE TABLE entrez_to_ensembl AS
SELECT
	CASE
		WHEN [NCBI.EntrezID] = '' THEN [HGNC.EntrezID]
		ELSE [NCBI.EntrezID]
	END AS EntrezID,
	CASE
		WHEN [Ensembl.EnsemblID] = '' THEN [HGNC.EnsemblID]
		ELSE [Ensembl.EnsemblID]
	END AS EnsemblID
FROM hgnc_raw
/* this looses 19060 - 18895 genes */
/* TODO find out which ones, and how I can avoid this */
WHERE EntrezID != '' AND EnsemblID != ''

