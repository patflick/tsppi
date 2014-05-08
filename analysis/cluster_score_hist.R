
library(ggplot2)
source("expr_utils.R")
source("plot_tiles.R")

get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all")
    return(exprs)
}

get_cluster_data <- function(ppi_name="string", expr_name="gene_atlas")
{
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')
    query <- paste("SELECT * FROM clustering_scoring_results ",
                   "WHERE size >= 4 AND ppi='", ppi_name,
                   "' AND expr='",expr_name, "'", sep="")
    data <- dbGetQuery(con, query)
    return(data)
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
get_top <- function(ppi_name="string", expr_name="gene_atlas")
{
    # TODO: top 10% PER PxExCxT
    data <- get_cluster_data(ppi_name, expr_name)
    data$mod_by_size <- data$modularity/data$size

    # per T
    df <- data.frame()
    types <- unique(data$type)
    for (t in types)
    {
        data_t <- data[which(data$type == t),]
        top10_thres <- quantile(data_t$modularity,0.8)
        data_top10 <- data_t[which(data_t$modularity >= top10_thres),]
        data_other <- data_t[which(data_t$modularity < top10_thres),]
        n_clusters <- dim(data_t)[1]
        mean_cl_size <- mean(data_t$size)

        tt <- t.test(data_top10$bpscore, data_other$bpscore)
        mean_top10 <- tt$estimate[1]
        mean_other <- tt$estimate[2]
        greater <- mean_top10 > mean_other
        t_pvalue <- tt$p.value
        ndf <- data.frame(ppi=ppi_name, expr=expr_name,
                          type=t,
                          n_clusters=n_clusters, mean_cl_size=mean_cl_size,
                          mean_top=mean_top10, mean_other=mean_other,
                          greater, pval=t_pvalue,
                          stringsAsFactors=FALSE)
        if (dim(df)[1] == 0) {
            df <- ndf
        } else {
            df <- rbind(df, ndf)
        }
    }
    return(df)
}

get_top_summary <- function(ppi_name="string", expr_name="gene_atlas")
{

    top_data <- get_top(ppi_name, expr_name)
    # can't do accumulated z-scores, because we don't have any
}

get_all_top <- function()
{
    df <- concat_ppi_expr_dataframe(get_top)
    return(df)
}
