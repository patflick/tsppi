

/* create table that has for each gene and tissue/cell.type
 * whether it is expressed or not (0: not expressed, 1: expressed,
 * no entry: no reliable data)
 */
DROP TABLE IF EXISTS hpa_tissue_expr_all;

CREATE TABLE hpa_tissue_expr_all AS
SELECT
	Gene,
	Tissue,
	[Cell.Type],
	CASE WHEN
		(
		Reliability = "High"
        OR
        Reliability = "Medium"
        OR
        Reliability = "Supportive"
        )
        AND
        (
	        NOT
	        (
	        Level = "Negative"
	        OR
	        Level = "None"
	        )
        )
        THEN 1
        ELSE 0
    END as Expressed
    FROM hpa_normal_tissue
    WHERE
	(
		Reliability = "High"
        OR
        Reliability = "Medium"
        OR
        Reliability = "Supportive"
    )
    ;
    

/* count each protein as expressed (==1) if
 * it is expressed in at least one cell type per tissue!
 */
DROP TABLE IF EXISTS hpa_tissue_expr;

CREATE TABLE hpa_tissue_expr AS
SELECT
	Gene,
	Tissue,
	max(Expressed) as Expressed
FROM  hpa_tissue_expr_all
GROUP BY Gene, Tissue

;


/* create tissue/cell.type specific expression
 * counts for edges in the PPI
 */

DROP TABLE IF EXISTS ppi_edge_expr;
CREATE TABLE ppi_edge_expr AS
SELECT
	a.Gene1,
	a.Gene2,
	SUM(CASE
		WHEN
			b.Expressed = 1 
			AND c.Expressed = 1
		THEN 1 
		ELSE 0 
		END) AS ExpressedCount,
	count() AS TotalCount
FROM ppi_genes AS a
INNER JOIN hpa_tissue_expr AS b
	ON a.Gene1 = b.Gene
INNER JOIN hpa_tissue_expr AS c
	ON a.Gene2 = c.Gene
	AND b.Tissue = c.Tissue
	/*AND b.[Cell.Type] = c.[Cell.Type]*/
GROUP BY
	a.Gene1, a.Gene2;

ALTER TABLE ppi_edge_expr ADD ExpressedFraction float;

UPDATE ppi_edge_expr
SET ExpressedFraction = ExpressedCount*1.0/TotalCount;





