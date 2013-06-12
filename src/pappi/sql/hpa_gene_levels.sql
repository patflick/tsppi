
DROP TABLE IF EXISTS hpa_gene_levels;
DROP TABLE IF EXISTS tmp_hpa_gene_levels;

CREATE TABLE tmp_hpa_gene_levels AS
SELECT
		a.Gene as Gene,
		count(CASE WHEN a.Level='High' OR a.Level="Strong" THEN 1 ELSE NULL END) as CountHigh,
		count(CASE WHEN a.Level='Medium' OR a.Level="Moderate" THEN 1 ELSE NULL END) as CountMedium,
		count(CASE WHEN a.Level='Low' OR a.Level="Weak" THEN 1 ELSE NULL END) as CountLow,
		count() as CountTotal
	FROM hpa_normal_tissue AS a
	/*
WHERE
	(
    a.[Expression.type] = "APE"
    AND
    (
        a.Reliability = "High"
        OR
        a.Reliability = "Medium" 
    )
    )
    
    OR
    (
    a.[Expression.type] = "Staining"
    AND
	a.Reliability = "Supportive"
    )
    
	*/
GROUP BY Gene
;

CREATE TABLE hpa_gene_levels AS
SELECT
	Gene,
	CountHigh,
	CountMedium,
	CountLow,
	CountHigh+CountMedium+CountLow as CountExpressed,
	CountTotal
FROM tmp_hpa_gene_levels
ORDER BY CountHigh+CountMedium+CountLow ASC;

DROP TABLE tmp_hpa_gene_levels

