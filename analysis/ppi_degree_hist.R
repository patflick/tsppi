# TODO: Add comment
# 
# Author: flick
###############################################################################
DATABASE_FILE = 'hpaDB.sqlite'

if (.Platform$OS.type == "unix")
{
	setwd("/cygdrive/d/PPI/");
} else {
	setwd("D:\\PPI");
}

dbpath = file.path(getwd(), DATABASE_FILE)

library('RSQLite')

drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname = dbpath)

data <- dbGetQuery(con, "
SELECT
		a.Gene,
		count(CASE WHEN a.Level='High' THEN 1 ELSE NULL END) as CountHigh,
		count(CASE WHEN a.Level='Medium' THEN 1 ELSE NULL END) as CountMedium,
		count(CASE WHEN a.Level='Low' THEN 1 ELSE NULL END) as CountLow,
		b.Degree as Degree
	FROM hpa_normal_tissue AS a
	INNER JOIN 
		ppi_node_degrees AS b
	ON a.Gene = b.[a.Gene]
WHERE
    [Expression.type] = 'APE'
	
    AND
    (
        Reliability = 'High'
        OR
        Reliability = 'Medium' 
    )
	
GROUP BY Gene
ORDER BY CountHigh+CountMedium+CountLow ASC
")


dbDisconnect(con)


# calculate the power law
library(igraph)
powerlaw <- power.law.fit(data$Degree)
print(powerlaw)



# save as png
#png("degree_distr.png", height=1050, width=1680)

par(mfrow=c(2,2))

thresholds<-c(0.2,0.4,0.6,0.8)

for (threshold in thresholds) {
	n <- nrow(data)
	degrees_high = data$Degree[1:floor(threshold*n)]
	degrees_low = data$Degree[(floor(threshold*n)+1):n]
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