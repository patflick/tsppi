# plots the histogram of protein expression in style of the plot
# for edge expression
# Author: flick
###############################################################################

# load utils
source("utils.R", chdir=TRUE)

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


# load ggplot2
library(ggplot2)

# choose bin width
binwidth <- 0.02

fig = ggplot(xlim=c(0,1))
fig = fig + stat_bin(data=expr_data, binwidth=binwidth, alpha=0.5, aes(x=ExpressedFraction, y=5*..count../sum(..count..), fill="Actual Data"), position="identity")
fig = fig + ylab("Fraction")
fig = fig + xlab("Fraction of cell types the protein is expressed")
fig = fig + labs(title="Histogram of Protein Expression")
fig = fig + scale_fill_manual(name="Protein Expression",values=c("blue"))
fig = fig + theme(legend.position=c(1,1),legend.justification=c(1,1))
#actually plot:
fig







