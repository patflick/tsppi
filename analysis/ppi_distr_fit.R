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


old_rand_subset_property <- function(ppi_graph, genes)
{
	n <- length(genes)
	# creates a random permutation of the genes
	threshold <- 0.5
	genes_subset <- sample(genes, floor(threshold*n))
	
	# get the induced subgraphs
	sub_ppi <- subset_ppi(ppi_graph, genes_subset)
	
	# calculate property of random subset
	#ass_total <- assortativity.degree(full_ppi_graph, directed=FALSE)
	ass_subset <- assortativity.degree(sub_ppi, directed=FALSE)
	return(ass_subset)
}

full_ppi_graph <- simplify(full_ppi_graph)

#vertex_property <- degree(full_ppi_graph)
#vertex_property <- betweenness(full_ppi_graph, directed=FALSE)
sh_paths = shortest.paths(full_ppi_graph)
sh_paths[is.infinite(sh_paths)] = NA
#vertex_property <- rowMeans(sh_paths, na.rm=TRUE)
graph_knn <- graph.knn(full_ppi_graph)
vertex_property <- graph_knn$knn


# get a random samples for assortativity
ass_stat <- c()
for (i in 1:10000){
    ass_stat <- c(ass_stat, rand_subset_property(full_ppi_graph, data$Gene, function(graph) assortativity(graph, types1=vertex_property[which(names(vertex_property) %in% V(graph)$name)], directed=FALSE)))
    if (i %% 2000 == 0){
        print(i)
    }
}


#for (threshold in thresholds) {
    
    threshold <- 0.5
    n <- nrow(data)
    # get the high and low specificity genes
    high_specificity_genes <- data$Gene[1:floor(threshold*n)]
    low_specificity_genes <- data$Gene[(floor(threshold*n)+1):n]
    
    # get the induced subgraphs
    low_spec_ppi <- subset_ppi(full_ppi_graph, low_specificity_genes)
    high_spec_ppi <- subset_ppi(full_ppi_graph, high_specificity_genes)
    
    # get node degrees of subgraphs
    #degrees_high = degree(high_spec_ppi)
    #degrees_low = degree(low_spec_ppi)
    #list[degrees_high, degrees_low] <- subgraph_vertex_property(full_ppi_graph, high_specificity_genes, low_specificity_genes, property_function=function(graph) betweenness(graph, directed=FALSE), method="global")
    
    ass <- assortativity.degree(full_ppi_graph, directed=FALSE)
    ass_low <- assortativity(low_spec_ppi, types1=vertex_property[which(names(vertex_property) %in% V(low_spec_ppi)$name)], directed=FALSE)
    ass_high <- assortativity(high_spec_ppi, types1=vertex_property[which(names(vertex_property) %in% V(high_spec_ppi)$name)], directed=FALSE)
    print(list(global=ass,low=ass_low,high=ass_high))
   

plot(density(ass_stat), main="Assortativity of node degree for random permutations")
abline(v = ass_low, col="blue")
abline(v = ass_high, col="red")



    #degrees_high=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh > threshold)]
    #degrees_low=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh <= threshold)]
#    plot(density(degrees_low), col="blue", log="x", main=paste(c("Degree distribution for threshold ", threshold),sep=""), xlab="Node degree")
#    lines(density(degrees_high),col="red")

plot(density(vertex_property[which(names(vertex_property) %in% V(high_spec_ppi)$name)]), col="red")#, log="xy")
lines(density(vertex_property[which(names(vertex_property) %in% V(low_spec_ppi)$name)]), col="blue")

low_spec_properties = vertex_property[which(names(vertex_property) %in% V(low_spec_ppi)$name)]
high_spec_properties = vertex_property[which(names(vertex_property) %in% V(high_spec_ppi)$name)]
#par(mfrow=c(2,1))
#hist(high_spec_properties, breaks=length(high_spec_properties))
#hist(low_spec_properties, breaks=length(low_spec_properties))

legend("topright",legend=c("High Specificity", "Low Specificity"), fill=c("red","blue"))
#}





dbDisconnect(con)

# close png output
#dev.off()
