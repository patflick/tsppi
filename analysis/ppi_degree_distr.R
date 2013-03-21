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
")






# save as png
png("degree_distr.png", height=800, width=1024)

par(mfrow=c(2,2))

thresholds<-c(15,30,45,60)

for (threshold in thresholds) {
	degrees_high=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh > threshold)]
	degrees_low=data$Degree[which(data$CountLow+data$CountMedium+data$CountHigh <= threshold)]
	plot(density(degrees_high), col="blue", main=paste(c("Degree distribution for threshold ", threshold),sep=""), xlab="Node degree")
	lines(density(degrees_low),col="red")
	legend("topright",legend=c("High Specificity", "Low Specificity"), fill=c("red","blue"))
}




dbDisconnect(con)

# close png output
dev.off()