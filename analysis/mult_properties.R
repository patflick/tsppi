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

thresholds <- c(0.3, 0.4)


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

# define number of graphs
par(mfcol=c(2,2))

for (threshold in thresholds)
{
    # get the high and low specificity genes
    list[low_specificity_genes, high_specificity_genes] <- get_genes_in_specificity_classes(threshold)


    source("plots/property_distribution.R", chdir=TRUE)

    # plot degree distr
    # vertex_property <- degree(full_ppi_graph)
    # get specific properties
    # list[low_spec_properties, high_spec_properties] <- get_vertex_properties_in_spec_classes(vertex_property, threshold)
    # plot_log_hist_lm(low_spec_properties, high_spec_properties, "Degree")
    # title(main=paste("Degree distr. lm, threshold ", 100*threshold, " %", sep="", collapse=""))

    source("plots/property_by_degree.R")
    vertex_property <- betweenness(full_ppi_graph, directed=FALSE)
    list[low_spec_properties, high_spec_properties] <- get_vertex_properties_in_spec_classes(vertex_property, threshold)
    mean_log_lm_plot_by_degree(full_ppi_graph, low_spec_properties, high_spec_properties, "Average betweeness")
    #boxplot_by_degree(full_ppi_graph, low_spec_properties, high_spec_properties, "Average betweeness")

    title(main=paste("Avg. Betweeness, threshold ", 100*threshold, " %", sep="", collapse=""))
    graph_knn <- graph.knn(full_ppi_graph)
    vertex_property <- graph_knn$knn
    list[low_spec_properties, high_spec_properties] <- get_vertex_properties_in_spec_classes(vertex_property, threshold)
    mean_log_lm_plot_by_degree(full_ppi_graph, low_spec_properties, high_spec_properties, "Average neighboor degree")

    title(main=paste("Avg. neighboor degree, threshold ", 100*threshold, " %", sep="", collapse=""))
    # sh_paths = shortest.paths(full_ppi_graph)
    # sh_paths[is.infinite(sh_paths)] = NA
    # vertex_property <- rowMeans(sh_paths, na.rm=TRUE)
    # list[low_spec_properties, high_spec_properties] <- get_vertex_properties_in_spec_classes(vertex_property, threshold)
    #  plot_log_hist(low_spec_properties, high_spec_properties, "Shortest path")

}


