# TODO: Add comment
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



# save as png
#png("degree_distr.png", height=800, width=1024)




full_ppi_graph <- simplify(full_ppi_graph)

#vertex_property <- degree(full_ppi_graph)
#vertex_property <- betweenness(full_ppi_graph, directed=FALSE)
#sh_paths = shortest.paths(full_ppi_graph)
#sh_paths[is.infinite(sh_paths)] = NA
#vertex_property <- rowMeans(sh_paths, na.rm=TRUE)
graph_knn <- graph.knn(full_ppi_graph)
#degrees <- degree(full_ppi_graph)

vertex_property <- graph_knn$knn[!is.na(graph_knn$knn)]
#non_nan_degrees <- degrees[match(names(non_nan_knn), names(degrees))]

# TODO average per degree


#deg_knn <- as.data.frame(cbind(non_nan_degrees, non_nan_knn))

#knn_deg_means <- aggregate(deg_knn["non_nan_knn"], by = deg_knn["non_nan_degrees"], FUN=mean)

#plot(non_nan_degrees, non_nan_knn, log="xy")


# TODO stddev of ALL data, meaning for every vertex, the degrees of all of his neighboors into a table
# TODO only then aggregate with mean and std-dev
# TODO especially no multi-stage mean (messes up the total std-dev)

#for (threshold in thresholds) {




threshold <- 0.5
# get the high and low specificity genes
high_specificity_genes <- data$Gene[1:floor(threshold*n)]
low_specificity_genes <- data$Gene[(floor(threshold*n)+1):n]

# get promiscuous and specific genes
source("gene_spec_classification.R", chdir=TRUE)
list[promiscuous_properties, specific_properties] <- get_vertex_properties_in_spec_classes(vertex_property, threshold)

# load plot function for plotting by degree
source("plots/property_by_degree.R", chdir=TRUE)

par(mfrow=c(3,1))

scatterplot_by_degree(full_ppi_graph, promiscuous_properties, specific_properties, "Average neighboor degree")
boxplot_by_degree(full_ppi_graph, promiscuous_properties, specific_properties, "Average neighboor degree")
mean_log_lm_plot_by_degree(full_ppi_graph, promiscuous_properties, specific_properties, "Average neigboor degree")



#degrees_high=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh > threshold)]
#degrees_low=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh <= threshold)]
#    plot(density(degrees_low), col="blue", log="x", main=paste(c("Degree distribution for threshold ", threshold),sep=""), xlab="Node degree")
#    lines(density(degrees_high),col="red")

#plot(density(vertex_property[which(names(vertex_property) %in% V(high_spec_ppi)$name)]), col="red")#, log="xy")
#lines(density(vertex_property[which(names(vertex_property) %in% V(low_spec_ppi)$name)]), col="blue")

#low_spec_properties = vertex_property[which(names(vertex_property) %in% V(low_spec_ppi)$name)]
#high_spec_properties = vertex_property[which(names(vertex_property) %in% V(high_spec_ppi)$name)]
#par(mfrow=c(2,1))
#hist(high_spec_properties, breaks=length(high_spec_properties))
#hist(low_spec_properties, breaks=length(low_spec_properties))

#legend("topright",legend=c("High Specificity", "Low Specificity"), fill=c("red","blue"))
#}





dbDisconnect(con)

# close png output
#dev.off()
