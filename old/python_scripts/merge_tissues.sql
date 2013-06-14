CREATE TABLE "merged_tissues" AS
SELECT
	a.Gene,
	a.Level AS Level1,
	c.Score AS Score1,
	b.Level AS Level2,
	d.Score AS Score2,
	c.Score - d.Score AS Diff
FROM tissue1 AS a 
INNER JOIN tissue2 AS b
	ON a.Gene = b.Gene
INNER JOIN LevelMapping AS c
	ON a.Level = c.Level
INNER JOIN LevelMapping AS d
	ON b.Level = d.Level
WHERE
	a."Expression.type" = "APE"
AND
	b."Expression.type" = "APE"
AND
(
	a.Reliability = "High"
	OR
	a.Reliability = "Medium"
)
AND
(
	b.Reliability = "High"
	OR
	b.Reliability = "Medium"
)