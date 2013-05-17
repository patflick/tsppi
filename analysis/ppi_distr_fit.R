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
# choose vertex property to calculate statistics about
########################################################

# 1.) vertex degree
 vertex_property <- degree(full_ppi_graph)

# 2.) betweenness centrality
# vertex_property <- betweenness(full_ppi_graph, directed=FALSE)

# 3.) avg shortest path
#sh_paths = shortest.paths(full_ppi_graph)
#sh_paths[is.infinite(sh_paths)] = NA
#vertex_property <- rowMeans(sh_paths, na.rm=TRUE)

# 4.) average neighboor degree
# graph_knn <- graph.knn(full_ppi_graph)
# vertex_property <- graph_knn$knn



########################################################
# choose thresholds
########################################################

thresholds <- c(0.2, 0.3, 0.5)


########################################################
# Code starts here
########################################################


# load the specific and promiscuous genes
source("gene_spec_classification.R", chdir=TRUE)


# get all genes for the random sampling
all_genes <- get_all_hpa_genes()


# define a function for the assortativity of subsets of the full graph
subset_assortativity <- function(full_graph, all_vertex_properties, vertex_subset)
{
    sub_graph <- subset_ppi(full_graph, vertex_subset)
  
    subset_properties <- all_vertex_properties[match(V(sub_graph)$name, names(all_vertex_properties))] 

    ass <- assortativity(sub_graph, types1=subset_properties, directed=FALSE)
    return(ass)
}

subset_internal_edge_percentage <- function(full_graph, subgraph_vertices, vertex_subsubset)
{
    graph <- subset_ppi(full_graph, subgraph_vertices)
    
    # filter the subsubset to only include vertices which are in the graph
    ss <- vertex_subsubset[which(vertex_subsubset %in% V(graph)$name)]

    # get all neighboors of the vertex_subsubset
    neighboors <- V(graph)[nei(ss)]
    neighboor_genes <- neighboors$name

    # check how many neighboors are again in the subsub set
    edge_internal <- neighboor_genes %in% ss
    percentage_internal_edges <- sum(edge_internal)/length(edge_internal)
    
    return(percentage_internal_edges)
}

# define number of graphs
par(mfcol=c(2,length(thresholds)))

for (threshold in thresholds)
{
    # get the high and low specificity genes
    list[low_specificity_genes, high_specificity_genes] <- get_genes_in_specificity_classes(threshold)

    # load the random permutation test plots
    source("plots/rand_perm_density.R", chdir=TRUE)
    # do random permutation and plot it according to the assortativity subset function
    #plot_rand_perm_density_2class(high_specificity_genes, low_specificity_genes, function(gene_subset) subset_assortativity(full_ppi_graph, vertex_property, gene_subset), num_perm_tests)

    # plot distribution for edge-internal-ness (meaning percentage of edges going
    # from promiscuous to promiscuous, rather then promiscuous->specific)
    # and the other way around
    plot_rand_perm_density_2class(
                high_specificity_genes,
                low_specificity_genes,
                function(gene_subset) subset_internal_edge_percentage(full_ppi_graph, all_genes, gene_subset),
                num_perm_tests, "Edge internalness")

    title(main=paste("Threshold ", threshold, " %", sep="", collapse=""))

    source("plots/property_distribution.R", chdir=TRUE)
    #plot_density(high_spec_properties, low_spec_properties, "Degree", TRUE)
    plot_log_hist_lm(high_spec_properties, low_spec_properties, "Degree")


}


