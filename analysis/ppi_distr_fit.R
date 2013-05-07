# does random permutation tests and calculates assortativity
# of various vertex properties
# 
# Author: flick
###############################################################################

# load utils
source("utils.R", chdir=TRUE)

# load sql config and get connection
source("sql_config.R", chdir=TRUE)
con <- get_sql_conn()

data <- dbGetQuery(con, "
SELECT
    Gene,
    CountHigh,
    CountMedium,
    CountLow,
    CountExpressed,
    CountTotal
FROM hpa_gene_levels
ORDER BY CountHigh+CountMedium+CountLow ASC

")



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
# vertex_property <- degree(full_ppi_graph)

# 2.) betweenness centrality
# vertex_property <- betweenness(full_ppi_graph, directed=FALSE)

# 3.) avg shortest path
#sh_paths = shortest.paths(full_ppi_graph)
#sh_paths[is.infinite(sh_paths)] = NA
#vertex_property <- rowMeans(sh_paths, na.rm=TRUE)

# 4.) average neighboor degree
graph_knn <- graph.knn(full_ppi_graph)
vertex_property <- graph_knn$knn



########################################################
# choose thresholds
########################################################

thresholds <- c(0.2, 0.3, 0.5)


par(mfcol=c(2,length(thresholds)))

for (threshold in thresholds) {
	
	
# get a random samples for assortativity
ass_stat_t_low <- c()
ass_stat_t_high <- c()
for (i in 1:num_perm_tests){
    ass_stat_t_low <- c(ass_stat_t_low, rand_subset_property(full_ppi_graph, data$Gene, function(graph) assortativity(graph, types1=vertex_property[which(names(vertex_property) %in% V(graph)$name)], directed=FALSE), subset_size=threshold))
	ass_stat_t_high <-  c(ass_stat_t_high, rand_subset_property(full_ppi_graph, data$Gene, function(graph) assortativity(graph, types1=vertex_property[which(names(vertex_property) %in% V(graph)$name)], directed=FALSE), subset_size=1-threshold))
    if (i %% 2000 == 0){
        print(i)
    }
}




n <- nrow(data)
# get the high and low specificity genes
high_specificity_genes <- data$Gene[1:floor(threshold*n)]
low_specificity_genes <- data$Gene[(floor(threshold*n)+1):n]

# get the induced subgraphs
low_spec_ppi <- subset_ppi(full_ppi_graph, low_specificity_genes)
high_spec_ppi <- subset_ppi(full_ppi_graph, high_specificity_genes)

# get the vertex property subsets for the two classes
low_spec_properties = vertex_property[which(names(vertex_property) %in% V(low_spec_ppi)$name)]
high_spec_properties = vertex_property[which(names(vertex_property) %in% V(high_spec_ppi)$name)]

# get assortativity for the subgraphs given by the promiscuous and specific proteins
ass_low <- assortativity(low_spec_ppi, types1=low_spec_properties, directed=FALSE)
ass_high <- assortativity(high_spec_ppi, types1=high_spec_properties, directed=FALSE)

den_below <- density(na.omit(ass_stat_t_low))
den_above <- density(na.omit(ass_stat_t_high))

xlims <- c(min(den_below$x,den_above$x, ass_low, ass_high), max(den_below$x,den_above$x, ass_low, ass_high))
ylims <- c(min(den_below$y,den_above$y), max(den_below$y,den_above$y))


# plot the density of the assortativity of the random sampling
plot(NaN, xlim=xlims, ylim=ylims, main="Assortativity", xlab="Assortativity", ylab="density of random samples")
lines(den_below, col="red")
lines(den_above, col="blue")
abline(v = ass_low, col="blue")
abline(v = ass_high, col="red")



    #degrees_high=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh > threshold)]
    #degrees_low=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh <= threshold)]
#    plot(density(degrees_low), col="blue", log="x", main=paste(c("Degree distribution for threshold ", threshold),sep=""), xlab="Node degree")
#    lines(density(degrees_high),col="red")






plot(density(na.omit(high_spec_properties)), main="Density", col="red")#, log="x")
lines(density(na.omit(low_spec_properties)), col="blue")

#par(mfrow=c(2,1))
#hist(high_spec_properties, breaks=length(high_spec_properties))
#hist(low_spec_properties, breaks=length(low_spec_properties))

legend("topright",legend=c("High Specificity", "Low Specificity"), fill=c("red","blue"))


}





dbDisconnect(con)

# close png output
#dev.off()
