# Loads the PPI from the database
# 
# Author: Patrick Flick
###############################################################################

load_ppi <- function() {
	# load sql config and get connection
	source("sql_config.R")
	con <- get_sql_conn()
	
	# load igraph library
	library(igraph)
	
	# load ppi network from db
	ppi_data <- dbGetQuery(con, "
	SELECT
		Gene1, Gene2
	FROM
		ppi_genes
	")

	# disconnect from db
	dbDisconnect(con)
	
	# load graph from edgelist
	ppi_graph <- graph.edgelist(as.matrix(ppi_data), directed=FALSE)
	
	# return the graph
	return(ppi_graph)
}

subset_ppi <- function(ppi_graph, genes)
{
	vertexes <- V(ppi_graph)
	ppi_subgraph <- induced.subgraph(ppi_graph, vertexes[which(vertexes$name %in% genes)], impl="create_from_scratch")
	return(ppi_subgraph)
}
