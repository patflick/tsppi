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



# save as png
#png("degree_distr.png", height=800, width=1024)

par(mfrow=c(2,2))

thresholds<-c(15,30,45,60)

for (threshold in thresholds) {
	# get the high and low specificity genes
	low_specificity_genes <- data$Gene[which(data$CountLow+data$CountMedium+data$CountHigh > threshold)]
	high_specificity_genes <- data$Gene[which(data$CountLow+data$CountMedium+data$CountHigh <= threshold)]
	
	# get the induced subgraphs
	low_spec_ppi <- subset_ppi(full_ppi_graph, low_specificity_genes)
	high_spec_ppi <- subset_ppi(full_ppi_graph, high_specificity_genes)
	
	# get node degrees of subgraphs
	degrees_high = degree(high_spec_ppi)
	degrees_low = degree(low_spec_ppi)
	
	#degrees_high=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh > threshold)]
	#degrees_low=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh <= threshold)]
	plot(density(degrees_low), col="blue", main=paste(c("Degree distribution for threshold ", threshold),sep=""), xlab="Node degree")
	lines(density(degrees_high),col="red")
	legend("topright",legend=c("High Specificity", "Low Specificity"), fill=c("red","blue"))
}




dbDisconnect(con)

# close png output
#dev.off()