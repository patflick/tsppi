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
		Gene,
		count(CASE WHEN Level='High' THEN 1 ELSE NULL END) as CountHigh,
		count(CASE WHEN Level='Medium' THEN 1 ELSE NULL END) as CountMedium,
		count(CASE WHEN Level='Low' THEN 1 ELSE NULL END) as CountLow
	FROM hpa_normal_tissue
WHERE
    [Expression.type] = 'APE'
	
    AND
    (
        Reliability = 'High'
        OR
        Reliability = 'Medium' 
    )
	
GROUP BY Gene
ORDER BY CountHigh+CountMedium+CountLow ASC, CountHigh ASC

")

dbDisconnect(con)


# save as png
png("hpa_tissue specificity.png", height=800, width=800)

# colors via http://colorschemedesigner.com/#3M11Tw0w0w0w0
colors=rev(c("#06246F", "#2040D0", "#5080FF"))
#colors=rev(topo.colors(3)[1:3])
plot(NA,ylim=c(0,max(data$CountHigh+data$CountMedium + data$CountLow)),xlim=c(0,dim(data)[1]), xlab="Genes", ylab="Cell Types")
# plot polygons as filled lines
legend("topleft", title="Staining Level", legend=c("Low","Medium","High"), fill=colors )
polygon(c(data$CountHigh+data$CountMedium + data$CountLow,0), col=colors[1], border=NA)
polygon(c(data$CountHigh+data$CountMedium,0), col=colors[2], border=NA)
polygon(c(data$CountHigh,0), col=colors[3], border=NA)

# close png output
dev.off()