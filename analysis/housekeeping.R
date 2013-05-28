



# load sql config and get connection
source("sql_config.R", chdir=TRUE)
con <- get_sql_conn()

data <- dbGetQuery(con, "
SELECT
    a.Gene as Gene,
    CountHigh,
    CountMedium,
    CountLow,
    CountExpressed,
    CountTotal,
    CASE WHEN b.Gene IS NULL THEN 0 ELSE 1 END AS HK_Gene,
    CASE WHEN c.Gene IS NULL THEN 0 ELSE 1 END AS TS_Gene
    
    FROM hpa_gene_levels as a
    LEFT OUTER JOIN hk_ensembl AS b
    ON a.Gene = b.Gene
    LEFT OUTER JOIN ts_ensembl AS c
    ON a.Gene = c.Gene
    ORDER BY CountHigh+CountMedium+CountLow ASC, CountHigh ASC
        ")


dbDisconnect(con)

colors=rev(c("#06246F", "#2040D0", "#5080FF"))
#colors=rev(topo.colors(3)[1:3])

x_range = c(0, dim(data)[1])
par(mar=c(5,4,4,5))
plot(NA,ylim=c(0,max(data$CountHigh+data$CountMedium + data$CountLow)),xlim=x_range, xlab="Gene Index", ylab="Cell Types")
# plot polygons as filled lines
polygon(c(data$CountHigh+data$CountMedium + data$CountLow,0), col=colors[1], border=NA)
polygon(c(data$CountHigh+data$CountMedium,0), col=colors[2], border=NA)
polygon(c(data$CountHigh,0), col=colors[3], border=NA)



# prepare the data for the chang et al. gene classification

hk_data <- density(which(data$HK_Gene == 1))
ts_data <- density(which(data$TS_Gene == 1))

hk_cum <- cumsum(data$HK_Gene)/max(cumsum(data$HK_Gene))
ts_cum <- cumsum(data$TS_Gene)/max(cumsum(data$TS_Gene))


par(new=TRUE)
plot(ts_data, col="red", type="l", lty=2, lwd=2, xaxt="n", xlim=x_range, yaxt="n", xlab="", ylab="", main="")
lines(hk_data, col="orange", type="l", lty=5, lwd=2)

# draw 50% lines
# abline(v=which(ts_cum > 0.5)[1], col="black", type="l", lty=3, lwd=2)
# abline(v=which(hk_cum > 0.5)[1], col="black", type="l", lty=3, lwd=2)

# get intersection between density plots:
intersX <- hk_data$x[as.logical(abs(diff(hk_data$y < ts_data$y)))]
 text(intersX, 70, "test", pos=4)
 abline(v=intersX, col="black", lty=3, lwd=2)

axis(4)
mtext("Chang et.al. classification density", side=4, line=3)

legend("bottomleft", title="Staining Level", legend=c("Low","Medium","High"), fill=colors )

legend("topright", legend=c("Tissue specific", "Housekeeping"), fill=c("red", "orange"))
title(main="Chang et. al. gene classification")


