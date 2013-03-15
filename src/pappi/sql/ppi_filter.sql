
/* 
 * filter out the backwards edges, 
 * making this an undirected instead of a directed graph
 * ( the edges in forward and backward direction are always scored identical
 * so no information is lost here)
 * 
 * */

DROP TABLE IF EXISTS ppi_unidir_proteins;
CREATE TABLE ppi_unidir_proteins
    ("id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, "Protein1" varchar(16), "Protein2" varchar(16), "Score" int);
    
INSERT INTO ppi_unidir_proteins
    ("Protein1", "Protein2", "Score")
SELECT
    Protein1,
    Protein2,
    Score
FROM ppi_raw_proteins
WHERE Protein1 < Protein2;


/*
 * Clean the ENSG to ENSP mapping (delete empty proteins)
 */
DELETE FROM ensg_to_ensp
WHERE ENSP = "";


/* Map Proteins to Genes */


DROP TABLE IF EXISTS ppi_genes;

CREATE TABLE ppi_genes (
    "id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
    "Gene1" varchar(16),
    "Gene2" varchar(16));
    
INSERT INTO ppi_genes
    (Gene1, Gene2)
SELECT DISTINCT
    CASE
        WHEN b.ENSG < c.ENSG 
        THEN b.ENSG
        ELSE c.ENSG
    END AS Gene1,
    CASE
        WHEN b.ENSG >= c.ENSG 
        THEN b.ENSG
        ELSE c.ENSG
    END AS Gene2
FROM ppi_unidir_proteins AS a
INNER JOIN
    ensg_to_ensp AS b
    ON a.Protein1 = b.ENSP
INNER JOIN
    ensg_to_ensp AS c
    ON a.Protein2 = c.ENSP
;



