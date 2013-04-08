
DROP TABLE IF EXISTS hpa_gene_levels;
DROP TABLE IF EXISTS tmp_hpa_gene_levels;

CREATE TABLE tmp_hpa_gene_levels AS
SELECT
		a.Gene as Gene,
		count(CASE WHEN a.Level='High' THEN 1 ELSE NULL END) as CountHigh,
		count(CASE WHEN a.Level='Medium' THEN 1 ELSE NULL END) as CountMedium,
		count(CASE WHEN a.Level='Low' THEN 1 ELSE NULL END) as CountLow,
		count() as CountTotal
	FROM hpa_normal_tissue AS a
WHERE
    [Expression.type] = 'APE'
	
    AND
    (
        Reliability = 'High'
        OR
        Reliability = 'Medium'
    )
	
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

