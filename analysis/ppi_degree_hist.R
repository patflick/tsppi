# TODO: Add comment
# 
# Author: flick
###############################################################################

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


# calculate the power law
library(igraph)
powerlaw <- power.law.fit(data$Degree)
print(powerlaw)



# save as png
#png("degree_distr.png", height=1050, width=1680)

par(mfrow=c(2,2))

thresholds<-c(0.2,0.4,0.6,0.8)

# TODO TODO
#  - differentiate between degree distribution of subgraphs
#    and the degrees in the full graph (both intersected with all APE genes, 
#    and ALL CCSB genes => different results, maybe the degree distribution
#    of the bigger (all CCSB genes) network is more special, the degree
#    distribution of only the subgraphs are a bit shitty because the networks are
#    too small (eg 472 nodes and 320 edges)
#  - compare subgraphs with all the other network analysis with igraph
#  - use also Staining (not only APE)

for (threshold in thresholds) {
	n <- nrow(data)
	# get the high and low specificity genes
	high_specificity_genes <- data$Gene[1:floor(threshold*n)]
	low_specificity_genes <- data$Gene[(floor(threshold*n)+1):n]
	
	# get the induced subgraphs
	low_spec_ppi <- subset_ppi(full_ppi_graph, low_specificity_genes)
	high_spec_ppi <- subset_ppi(full_ppi_graph, high_specificity_genes)
	
	# get node degrees of subgraphs
	degrees_high = degree(high_spec_ppi)
	degrees_low = degree(low_spec_ppi)
	
	#degrees_high = data$Degree[1:floor(threshold*n)]
	#degrees_low = data$Degree[(floor(threshold*n)+1):n]
	hist_high = hist(degrees_high, plot=FALSE,breaks=max(degrees_high))
	hist_low = hist(degrees_low, plot=FALSE,breaks=max(degrees_low))
	y1 = c(0,hist_low$density,0)
	y2 = c(0,hist_high$density,0) 
	plot(c(0,hist_low$breaks), y1 ,type="s", col="blue", xlim= c(0,max(c(hist_low$breaks,hist_high$breaks))),ylim= c(0,max(c(y1,y2))),main=paste(c("Degree distribution for threshold ", threshold),sep=""), xlab="Node degree", ylab="Density")
	lines(x=c(0,hist_high$breaks), y = y2, type="s", col="red")
	legend("topright",legend=c("High Specificity", "Low Specificity"), fill=c("red","blue"))
	pl_high <- power.law.fit(degrees_high)
	pl_low <- power.law.fit(degrees_low)
	text(100,0.4, paste(c("High Specificity\nPower Law: ",pl_high$alpha, "\nLow Specificity\nPower Law: ", pl_low$alpha), sep="", collapse=""), pos=4)
}






# close png output
#dev.off()