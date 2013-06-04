# does random permutation tests and calculates assortativity
# of various vertex properties
# 
# Author: flick
###############################################################################


# load full PPI network
source("ppi_network_functions.R", chdir=TRUE)
full_ppi_graph <- load_ppi()
full_ppi_graph <- simplify(full_ppi_graph)


########################################################
# choose number of random samples (permutation test)
#  -> more = more accurate, less = faster
########################################################
num_perm_tests <- 10000


########################################################
# Code starts here
########################################################


# load the specific and promiscuous genes
source("gene_spec_classification.R", chdir=TRUE)


# get all genes for the random sampling
all_genes <- get_all_hpa_genes()

# returns the number of edges in the induced subgraph
subgraph_number_edges <- function(full_graph, subgraph_vertices)
{
    graph <- subset_ppi(full_graph, subgraph_vertices)
    return(length(E(graph)))
}

# returns the number of vertices in the induced subgraph
subgraph_number_vertices <- function(full_graph, subgraph_vertices)
{
	graph <- subset_ppi(full_graph, subgraph_vertices)
	return(length(V(graph)))
}

# define number of graphs
par(mfcol=c(1,1))

NUM_THRESHOLDS <- 20

thresholds <- 0:NUM_THRESHOLDS/(NUM_THRESHOLDS)

num_edges_low_spec <- c()
num_edges_high_spec <- c()
num_vertices_low <- c()
num_vertices_high <- c()

for (threshold in thresholds)
{
    # get the high and low specificity genes
    list[low_specificity_genes, high_specificity_genes] <- get_genes_in_specificity_classes(threshold)

	num_edges_low_spec <- c(num_edges_low_spec, subgraph_number_edges(full_ppi_graph, low_specificity_genes))
	num_edges_high_spec <- c(num_edges_high_spec, subgraph_number_edges(full_ppi_graph, high_specificity_genes))
	num_vertices_low <- c(num_vertices_low,  subgraph_number_vertices(full_ppi_graph, low_specificity_genes))
	num_vertices_high <- c(num_vertices_high,  subgraph_number_vertices(full_ppi_graph, high_specificity_genes))
}


library(ggplot2)


fig = ggplot(data=data.frame(
		x=thresholds,
		edges_low=num_edges_low_spec,
		edges_high=num_edges_high_spec,
		vert_low=num_vertices_low,
		vert_high=num_vertices_high)
	, aes(x=x))
fig = fig + geom_line(aes(y=edges_low, color="|E| in promiscuous subgraph"))
fig = fig + geom_line(aes(y=edges_high, color="|E| in specific subgraph"))
fig = fig + geom_line(aes(y=vert_low, color="|V| in promiscuous subgraph"))
fig = fig + geom_line(aes(y=vert_high, color="|V| in specific subgraph"))
fig


