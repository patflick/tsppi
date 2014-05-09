
library(ggplot2)
source("expr_utils.R")
source("plot_tiles.R")

get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all")
    return(exprs)
}

clusterers <- c('PLP','PLM-gamma-1.0', 'PLM-gamma-5.0', 'PLM-gamma-10.0', 'PLM-gamma-50.0', 'PLM-gamma-100.0')



get_cluster_data <- function(ppi_name="string", expr_name="gene_atlas", clusterer=clusterers[1])
{
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')
    query <- paste("SELECT * FROM clustering_scoring_results ",
                   "WHERE size >= 4 AND ppi='", ppi_name,
                   "' AND expr='",expr_name, "'",
                   " AND clusterer = '", clusterer, "'",sep="")
    data <- dbGetQuery(con, query)
    return(data)
}

get_clusterer_data <- function(ppi_name="string", expr_name="gene_atlas")
{
    type <- "Global"
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')
    query <- paste("SELECT ppi, expr, clusterer, SUM(modularity) as total_mod, COUNT(), MIN(size), AVG(size), MAX(size) FROM clustering_scoring_results ",
                   "WHERE",
                   " size >= 4",
                   #" AND ppi='", ppi_name,"'",
                   #" AND expr='",expr_name, "'",
                   " AND type = '", type, "'",
                   " GROUP BY ppi, expr, clusterer", sep="")
    data <- dbGetQuery(con, query)
    return(data)
}


plot_cluster_distr <- function(ppi_name="string", expr_name="gene_atlas")
{
    type <- "GLOBAL"
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')
    query <- paste("SELECT * FROM clustering_scoring_results ",
                   "WHERE ppi='", ppi_name,
                   "' AND expr='",expr_name, "'",
                   " AND type='", type, "'",
                   " ORDER BY clusterer, size", sep="")
    data <- dbGetQuery(con, query)
    data$size_discrete <- factor(as.character(1:dim(data)[1]))
    data$clusterer_fac <- factor(data$clusterer, levels=rev(clusterers))

    fig <- ggplot(data, aes(clusterer_fac, size, fill=size_discrete)) +
        geom_bar(stat="identity", position = "stack") +
        coord_flip() +
        scale_fill_manual(values = rep(c("grey20", "grey40", "grey60"), ceiling(dim(data)[1]/3)), guide=FALSE) +
        xlab("Clustering Algorithm") +
        ylab("Cluster sizes") +
        labs(title="Cluster sizes by clustering algorithm")
        #scale_fill_hue()

    return(fig)
}

save_cluster_sizes <- function(ppi_name="string", expr_name="gene_atlas")
{
    # global ts network TS betweenness
    pdf(paste("../figs/cluster_size_distr.pdf", sep=""), width=8, height=3.5)
    p <- plot_cluster_distr(ppi_name, expr_name)
    plot(p)
    dev.off()

}

# get the top10 sizes of clusters per clusterer
get_top_cluster_distr <- function(ppi_name="string", expr_name="gene_atlas")
{
    type <- "GLOBAL"
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')
    query <- paste("SELECT * FROM clustering_scoring_results ",
                   "WHERE ppi='", ppi_name,
                   "' AND expr='",expr_name, "'",
                   " AND type='", type, "'",
                   " ORDER BY clusterer, size", sep="")
    data <- dbGetQuery(con, query)

    df <- data.frame()

    for (c in clusterers)
    {
        df <- rbind(df, tail(data[which(data$clusterer == c),], 10))
    }

    print(df)

    for (c in clusterers)
    {
        data_c <- data[which(data$clusterer == c),]
        print(paste(c, " num of clusters < 4 = ", sum(data_c$size[which(data_c$size < 4)])))
    }

    return(df)
}

plot_scores_box <- function()
{
    data <- get_cluster_data()
    fig <- ggplot(data, aes(x=type, y=bpscore, group)) +
            geom_boxplot() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

get_corrs <- function()
{
    df <- data.frame(ppi=character(0), expr=character(0), corr=double(0), stringsAsFactors=FALSE)
    for (p in get_ppis())
    {
        for (e in get_exprs())
        {
            data <- get_cluster_data(p, e)
            c <- cor.test(data$modularity/data$size, data$bpscore)
            new_row <- data.frame(ppi=p, expr=e, corr=c$estimate,stringsAsFactors=FALSE)
            df <- rbind(df, new_row)
        }
    }
    return(df)
}

concat_ppi_expr_dataframe <- function(func, ...)
{
    # result:
    df <- NULL

    ppis <- get_ppis()
    exprs <- get_exprs()
    # for all ppis and expression data sets
    for (p in ppis)
    {
        for (e in exprs)
        {
            # print out?
            print(paste("getting output for: ",p, " and ", e))
            #call the given function
            result <- func(p,e, ...)

            if (is.null(df))
            {
                # first call
                df <- result
            }
            else
            {
                # concat
                df <- rbind(df, result)
            }
        }
    }
    return (df)
}

# get top 10% of clusters
get_top <- function(ppi_name="string", expr_name="gene_atlas", clusterer=clusterers[1], types=NA)
{
    # TODO: top 10% PER PxExCxT
    data <- get_cluster_data(ppi_name, expr_name, clusterer)
    data$mod_by_size <- data$modularity/data$size

    # per T
    df <- data.frame()
    if (is.na(types))
    {
        types <- unique(data$type)
    }
    for (t in types)
    {
        data_t <- data[which(data$type == t),]
        n_clusters <- dim(data_t)[1]
        mean_cl_size <- mean(data_t$size)
        # sort by columns
        data_t <- data_t[with(data_t, order(-modularity)), ]
        top20_idx <- ceiling(dim(data_t)[1] * 0.2)
        data_top20 <- data_t[1:top20_idx,]
        data_other <- data_t[(top20_idx+1):(dim(data_t)[1]),]
        if (dim(data_top20)[1] < 2 || dim(data_other)[1] < 2)
        {
            # too few for statistics
            ndf <- data.frame(ppi=ppi_name, expr=expr_name,
                              type=t,
                              n_clusters=n_clusters, mean_cl_size=mean_cl_size,
                              mean_top=NA, mean_other=NA,
                              greater=FALSE, pval=NA,
                              stringsAsFactors=FALSE)
        } else {
            # get statistics
            tt <- t.test(data_top20$bpscore, data_other$bpscore)
            mean_top10 <- tt$estimate[1]
            mean_other <- tt$estimate[2]
            greater <- mean_top10 > mean_other
            t_pvalue <- tt$p.value
            ndf <- data.frame(ppi=ppi_name, expr=expr_name,
                              type=t,
                              n_clusters=n_clusters, mean_cl_size=mean_cl_size,
                              mean_top=mean_top10, mean_other=mean_other,
                              greater=greater, pval=t_pvalue,
                              stringsAsFactors=FALSE)
        }
        if (dim(df)[1] == 0) {
            df <- ndf
        } else {
            df <- rbind(df, ndf)
        }
    }
    return(df)
}

get_top_global <- function()
{
    # use the PLM-gamma-50 clustering
    clusterer <- "PLM-gamma-50.0"
    expr <- "rnaseq_atlas"

    df <- data.frame()

    for (p in get_ppis())
    {
        top_data <- get_top(p, expr, clusterer, "EdgeScoring")
        df <- rbind(df, top_data)
    }

    # make the table look nice
    df <- df[,c("ppi", "n_clusters", "mean_cl_size", "mean_top", "mean_other", "greater", "pval")]
    df$ppi <- to_short_ppi_name(df$ppi)


    return(df)
}

get_all_top <- function(clusterer=clusterers[1])
{
    df <- concat_ppi_expr_dataframe(get_top,clusterer)
    return(df)
}
