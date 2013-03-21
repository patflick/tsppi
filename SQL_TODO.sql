
CREATE TABLE ppi_genes_bi AS SELECT Gene1, Gene2 FROM ppi_genes UNION SELECT Gene2 as Gene1, Gene1 as Gene2 FROM ppi_genes

DROP TABLE IF EXISTS ppi_node_degrees
CREATE TABLE ppi_node_degrees AS SELECT a.Gene, count() as Degree FROM (SELECT DISTINCT Gene FROM hpa_ape_reliable_score) AS a INNER JOIN ppi_genes_bi as b ON a.Gene=b.Gene1 GROUP BY a.Gene

SELECT sum(Degree) FROM ppi_node_degrees
