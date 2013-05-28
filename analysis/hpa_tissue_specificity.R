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
ORDER BY CountHigh+CountMedium+CountLow ASC, CountHigh ASC

")

dbDisconnect(con)


# save as png
#png("hpa_tissue specificity.png", height=800, width=800)

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
# dev.off()
