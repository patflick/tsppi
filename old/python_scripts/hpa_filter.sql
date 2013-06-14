/* Filters the HPA data */
DROP TABLE IF EXISTS hpa_ape_reliable_score;

/* create score mapping table */

DROP TABLE IF EXISTS hpa_ape_score_mapping;
CREATE TABLE hpa_ape_score_mapping ("Level" varchar(16), "Score" float);
INSERT INTO hpa_ape_score_mapping VALUES ("High", 3.0);
INSERT INTO hpa_ape_score_mapping VALUES ("Medium", 2.0);
INSERT INTO hpa_ape_score_mapping VALUES ("Low", 1.0);
INSERT INTO hpa_ape_score_mapping VALUES ("None", 0.0);


/* filter HPA data to only contain APE and High or Medium Reliability */
/* and add the numeric score */
CREATE TABLE hpa_ape_reliable_score
AS
SELECT 
    a.Gene,
    a.Tissue,
    a.[Cell.type],
    a.Level,
    a.[Expression.type],
    a.Reliability,
    b.Score
FROM hpa_normal_tissue AS a
INNER JOIN
    hpa_ape_score_mapping AS b
    ON a.Level = b.Level
WHERE
    a.[Expression.type] = "APE"
    AND
    (
        a.Reliability = "High"
        OR
        a.Reliability = "Medium" 
    )
;


/* create a summary table */

DROP TABLE IF EXISTS hpa_ape_reliable_tissues;
CREATE TABLE hpa_ape_reliable_tissues
	("id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, "Tissue" varchar(256), "Cell.Type" varchar(256), "NrProteins" int);

INSERT INTO hpa_ape_reliable_tissues
	(Tissue, [Cell.Type], NrProteins)
    SELECT Tissue,
           [Cell.type] as [Cell.type],
           count() AS NrProteins
        FROM hpa_ape_reliable_score
    GROUP BY Tissue,
              [Cell.type];

              
/* filter cell types for cutoff value */

DROP TABLE IF EXISTS selected_cell_types;
CREATE TABLE selected_cell_types
AS
    SELECT Tissue, [Cell.Type]
    FROM hpa_ape_reliable_tissues
    WHERE NrProteins >= 2080
    ;

    
/* create mean tissue expression table */

DROP TABLE IF EXISTS hpa_mean_tissue;

CREATE TABLE hpa_mean_tissue
          AS SELECT a.Gene as Gene,
                    avg(a.Score) AS Score,
                    count() AS Count
               FROM hpa_ape_reliable_score AS a
               INNER JOIN selected_cell_types AS b
                   ON a.Tissue = b.Tissue
                       AND
                      a.[Cell.type] = b.[Cell.type]
              GROUP BY Gene
              HAVING count() = (SELECT count() FROM selected_cell_types)
;


/* insert mean tissue and celltype into the reliable_score table */

INSERT INTO hpa_ape_reliable_score
	(Gene, Tissue, [Cell.type], Level, [Expression.type], Reliability, Score)
SELECT 
    Gene,
    "mean",
    "mean",
    "undefined",
    "undefined",
    "undefined",
    Score
FROM hpa_mean_tissue;


/* update summary table */

DELETE FROM hpa_ape_reliable_tissues;
INSERT INTO hpa_ape_reliable_tissues
	(Tissue, [Cell.Type], NrProteins)
    SELECT Tissue,
           [Cell.type] as [Cell.type],
           count() AS NrProteins
        FROM hpa_ape_reliable_score
    GROUP BY Tissue,
              [Cell.type];
