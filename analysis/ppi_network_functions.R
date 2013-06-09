# Loads the PPI from the database
# 
# Author: Patrick Flick
###############################################################################

# load igraph library
library(igraph)

load_ppi <- function(table_name="ppi_hgnc") {
	# load sql config and get connection
	source("sql_config.R")
	con <- get_sql_conn()
	
	# load ppi network from db
	ppi_data <- dbGetQuery(con, paste(c("
	SELECT
		Gene1, Gene2
	FROM ", table_name), sep="", collapse=""))

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
	ppi_subgraph <- induced.subgraph(ppi_graph, vertexes[which(vertexes$name %in% genes)], impl="auto")
	return(ppi_subgraph)
}


rand_subset_property <- function(ppi_graph, genes, property_function=function(graph) assortativity.degree(graph, directed=FALSE),subset_size = 0.5)
{
    # TODO wtf mache ich mit dem hier? wieso klappt es nicht?
#    ppi_merged <- subset_ppi(ppi_graph, genes)
 #   genes <- V(ppi_merged)


	n <- length(genes)
	# creates a random permutation of the genes
	genes_subset <- sample(genes, floor(subset_size*n))
	
	# get the induced subgraphs
	sub_ppi <- subset_ppi(ppi_graph, genes_subset)
	
	# calculate property of random subset
	#ass_total <- assortativity.degree(full_ppi_graph, directed=FALSE)
	subset_property <- property_function(sub_ppi)
	return(subset_property)
}


subgraph_vertex_property <- function(ppi_graph, x, y, property_function=degree, method="global")
{
	
	if (method == "global")
	{
		# the network property of the whole network is calculated
		# and then the results for the two sets returned
		vertex_properties_all <- property_function(ppi_graph)
		vertex_properties_x <- vertex_properties_all[intersect(x, names(vertex_properties_all))]
		vertex_properties_y <- vertex_properties_all[intersect(y, names(vertex_properties_all))]
	}
	else if (method == "subgraph")
	{
		# create subgraphs, then calculate the vertex properties on them
		subgraph_x <- subset_ppi(ppi_graph, x)
		vertex_properties_x <- property_function(subgraph_x)
		
		subgraph_y <- subset_ppi(ppi_graph, y)
		vertex_properties_y <- property_function(subgraph_y)
	}
	else if (method == "combined_subgraph")
	{
		# use the subgraph defined by combining x and y for vertex properties
		subgraph_xy <- subset_ppi(ppi_graph, c(x,y))
		vertex_properties_xy <- property_function(subgraph_xy)
		vertex_properties_x <- vertex_properties_xy[intersect(x, names(vertex_properties_xy))]
		vertex_properties_y <- vertex_properties_xy[intersect(y, names(vertex_properties_xy))]
	}
	else
	{
		stop("invalid argument method to `subgraph_vertex_property`")
	}
	
	return(list(x=vertex_properties_x, y=vertex_properties_y))
}
