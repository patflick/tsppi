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
		[Expression.type] || ' - ' || Reliability as Scoring,
		Tissue || ' / ' || [Cell.type] as [Tissue.Cell]
 	FROM hpa_normal_tissue
")

dbDisconnect(con)


# create box plot

tbl <- table(data$Scoring, data$Tissue.Cell)

j <- 1
for (col in colnames(tbl))
{
	i <- 1
	for (row in rownames(tbl))
	{
		x <- which(data$Tissue.Cell == col & data$Scoring == row)
		tbl[i,j] = data$Count[x]
		i <- i + 1
	}
	
	j <- j + 1
}

tbl <- tbl[c("APE - High", "APE - Medium", "APE - Low", "APE - Very low", "Staining - Supportive", "Staining - Uncertain"),]

colors <- c("#004400", "#008800", "#22cc22", "#77ff77", "#0000aa", "#5555ff")


png("hpa_data.png", height=1200, width=600)
oldmar <-par()$mar
oldmar[2] <- 25
par(mar=oldmar)
par(las=2)
barplot(tbl, horiz=TRUE, col=colors, xlab="Number of Genes")
par(xpd=TRUE)
#legend("topleft",inset=c(-0.4,0),legend=rownames(tbl), col=colors, fill=colors)
legend("bottomleft",inset=c(-0.4,0),legend=rownames(tbl), col=colors, fill=colors)
#axis(3)
dev.off()