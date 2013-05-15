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


# load gene classification
source("gene_spec_classification.R", chdir=TRUE)


par(mfrow=c(2,length(thresholds)))

for (threshold in thresholds) {

# get genes in their classes
list[low_spec_properties, high_spec_properties] <- get_vertex_properties_in_spec_classes(vertex_property, threshold)

h_nbreaks = 7

h_high <- hist(log10(high_spec_properties), breaks=h_nbreaks, plot=FALSE)
h_low <- hist(log10(low_spec_properties), breaks=h_nbreaks, plot=FALSE)

x_high <- 10**h_high$mids
y_high <- h_high$counts/sum(h_high$counts)

x_low <- 10**h_low$mids
y_low <- h_low$counts/sum(h_low$counts)



#plot(density(na.omit(high_spec_properties)), main="Density", col="red", log="xy")
#lines(density(na.omit(low_spec_properties)), col="blue")
plot(x_high, y_high, col="red", type="p", log="xy", xlim=c(1,max(x_high,x_low)), ylim=c(min(y_high, y_low),max(y_high,y_low)), pch=18)
abline(lm(log10(y_high)~log10(x_high)),col="red", lty=2)
lines(x_low, y_low, col="blue", type="p", pch=17)
abline(lm(log10(y_low)~log10(x_low)),col="blue",lty=2)


#par(mfrow=c(2,1))
#hist(high_spec_properties, breaks=length(high_spec_properties))
#hist(low_spec_properties, breaks=length(low_spec_properties))

legend("topright",legend=c("High Specificity", "Low Specificity"), fill=c("red","blue"))


}


# close png output
#dev.off()
