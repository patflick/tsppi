# plots the histogram of protein expression in style of the plot
# for edge expression
# Author: flick
###############################################################################

# load utils
source("utils.R", chdir=TRUE)

# load full PPI network
source("ppi_network_functions.R", chdir=TRUE)
ccsb_ppi_graph <- load_ppi("ppi_ccsb_hgnc")
mmc_ppi_graph <- load_ppi("ppi_mmc_hgnc")
#full_ppi_graph <- simplify(full_ppi_graph)


# load sql config and get connection
source("sql_config.R", chdir=TRUE)
con <- get_sql_conn()


expr_data <- dbGetQuery(con, "
	SELECT Gene, 
		SUM(Expressed) as CountExpressed,
		count() as CountTotal,
		SUM(Expressed)*1.0/count() as ExpressedFraction
	FROM hpa_tissue_expr
	GROUP BY Gene")


dbDisconnect(con)



ccsb_sub_graph <- subset_ppi(ccsb_ppi_graph, expr_data$Gene)
mmc_sub_graph <- subset_ppi(mmc_ppi_graph, expr_data$Gene)

# filter genes that are in the PPI
ccsb_expr_data <- expr_data[which(expr_data$Gene %in% V(ccsb_sub_graph)$name),]
mmc_expr_data <-  expr_data[which(expr_data$Gene %in% V(mmc_sub_graph)$name),]


# load ggplot2
library(ggplot2)

NUM_BREAKS <- 20

hist_breaks <- 0:NUM_BREAKS/NUM_BREAKS
h_expr <- hist(expr_data$ExpressedFraction, breaks=hist_breaks, plot=FALSE)
h_frac <- h_expr$counts/sum(h_expr$counts)
h_ccsb_expr <- hist(ccsb_expr_data$ExpressedFraction, breaks=hist_breaks, plot=FALSE)
h_ccsb_frac <- h_ccsb_expr$counts/sum(h_expr$counts)
h_mmc_expr <- hist(mmc_expr_data$ExpressedFraction, breaks=hist_breaks, plot=FALSE)
h_mmc_frac <- h_mmc_expr$counts/sum(h_expr$counts)


#
fig = ggplot(xlim=c(0,1))
fig = fig + geom_bar(data=data.frame(x=h_expr$mids, y=h_frac), alpha=0.5, aes(x=x,y=y, fill="All reliable HPA genes"),stat="identity")

fig = fig + geom_bar(data=data.frame(x=h_ccsb_expr$mids, y=h_ccsb_frac), aes(x=x, y=y, fill="PPI genes"), stat="identity")

#fig = fig + stat_bin(data=expr_data, alpha=0.5, aes(x=ExpressedFraction, y=5*..count../sum(..count..), fill="Actual Data"), position="identity")
fig = fig + ylab("Fraction")
fig = fig + xlab("Fraction of cell types the protein is expressed")
fig = fig + labs(title="Histogram of Protein Expression")
fig = fig + scale_fill_manual(name="Protein Expression",values=c("blue","orange"))
fig = fig + theme(legend.position=c(0,1),legend.justification=c(0,1))



fig2 = ggplot(xlim=c(0,1), data=data.frame(x=h_expr$mids, y1=h_ccsb_frac/h_frac, y2=h_mmc_frac/h_frac))
fig2 = fig2 + geom_point(aes(x=x, y=y1, color="CCSB"), shape=1)
fig2 = fig2 + geom_smooth(aes(x=x, y=y1, color="CCSB"), method=lm)
fig2 = fig2 + geom_point(aes(x=x, y=y2, color="MMC"), shape=2)
fig2 = fig2 + geom_smooth(aes(x=x, y=y2, color="MMC"), method=lm)



par(mfrow=c(1,2))
fig
fig2






