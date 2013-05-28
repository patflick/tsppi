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
num_perm_tests <- 1000


########################################################
# choose thresholds
########################################################

thresholds <- c(0.3)


########################################################
# Code starts here
########################################################


# load the specific and promiscuous genes
source("gene_spec_classification.R", chdir=TRUE)


# get all genes for the random sampling
all_genes <- get_all_hpa_genes()


# calc vertex degree
vertex_property <- degree(full_ppi_graph)




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

    # get specific properties
    list[low_spec_properties, high_spec_properties] <- get_vertex_properties_in_spec_classes(vertex_property, threshold)

    # load the random permutation test plots
    source("plots/rand_perm_density.R", chdir=TRUE)

    # do random permutation and plot it according to the assortativity subset function
    plot_rand_perm_density_2class(high_specificity_genes, low_specificity_genes, function(gene_subset) subset_assortativity(full_ppi_graph, vertex_property, gene_subset), num_perm_tests, "Assortativity (wrt. Degree)")


    title(main=paste("Assortativity, threshold ", 100*threshold, " %", sep="", collapse=""))

    source("plots/property_distribution.R", chdir=TRUE)
    
    plot_density(low_spec_properties, high_spec_properties, "Degree", FALSE)
    title(main=paste("Degree density, threshold ", 100*threshold, " %", sep="", collapse=""))

    plot_density(low_spec_properties, high_spec_properties, "Degree", TRUE)
    title(main=paste("Degree density (log), threshold ", 100*threshold, " %", sep="", collapse=""))

    plot_log_hist_lm(low_spec_properties, high_spec_properties, "Degree")
    title(main=paste("Degree distr. lm, threshold ", 100*threshold, " %", sep="", collapse=""))


}


