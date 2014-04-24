
library(ggplot2)
library(xtable) # for converting a data frame into a latex table
library(gridExtra) # for `grid.arrange`
library(igraph) # for power.law.fit

# load PPI name mapping
source("ppi_utils.R")

# returns node properties for ts graphs
get_node_properties <- function(ppi_name = "string", expr_name = "gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    ppi_node_prop_tbl <- paste(ppi_name,expr_name, "node_properties", sep="_")
    query <- paste("SELECT * FROM ", ppi_node_prop_tbl)
    node_props <- dbGetQuery(con, query)

    return(node_props)
}

get_graph_properties <- function(ppi_name="string", expr_name="gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    ppi_prop_tbl <- paste(ppi_name, expr_name, "properties", sep="_")
    query <- paste("SELECT * FROM ", ppi_prop_tbl)
    props <- dbGetQuery(con, query)

    return (props)
}

# get tissue specific node properties
get_ts_node_properties <- function(ppi_name = "string", expr_name = "gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    ppi_node_prop_tbl <- paste(ppi_name,expr_name, "ts_node_properties", sep="_")
    query <- paste("SELECT * FROM ", ppi_node_prop_tbl)
    node_props <- dbGetQuery(con, query)
    return(node_props)
}

# get node properties joined to expression counts (for HK, TS classification)
get_node_properties_exprcount <- function(ppi_name="string",expr_name="gene_atlas")
{
    # load sql config and get connection
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # create table names
    ppi_node_prop_tbl <- paste(ppi_name,expr_name, "node_properties", sep="_")
    expr_count_table <- paste(expr_name,"core_expr_counts", sep="_")

    # create SQL query to get ExpressionCount vs Degree
    query <- paste("SELECT a.*, b.ExpressedCount, b.TotalCount FROM ",
                   ppi_node_prop_tbl, " AS a INNER JOIN ", expr_count_table,
                   " AS b ON a.Gene = b.Gene WHERE ExpressedCount > 0")

    # load ppi network from db
    data <- dbGetQuery(con, query)
    return(data);
}

get_ts_node_properties_exprcount <- function(ppi_name="string",expr_name="gene_atlas")
{
    # load sql config and get connection
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # create table names
    ppi_node_prop_tbl <- paste(ppi_name,expr_name, "ts_node_properties", sep="_")
    expr_count_table <- paste(expr_name,"core_expr_counts", sep="_")

    # create SQL query to get ExpressionCount vs Degree
    query <- paste("SELECT a.*, b.ExpressedCount, b.TotalCount FROM ",
                   ppi_node_prop_tbl, " AS a INNER JOIN ", expr_count_table,
                   " AS b ON a.Gene = b.Gene WHERE ExpressedCount > 0")

    # load ppi network from db
    data <- dbGetQuery(con, query)
    return(data);
}

# plot the distribution of the betweenness centrality measures
plot_betweenness_distr <- function(ppi_name = "string", expr_name = "gene_atlas")
{
    data <- get_node_properties(ppi_name, expr_name)

    # get the betweenness frequency distribution
    bw <- as.integer(round(data$Betweenness, 100))
    bw_distr <- as.data.frame(prop.table(table(bw)))
    bw_distr$Betweenness <- as.integer(levels(bw_distr$bw))

    # plot on log-log scale (since it appears to be exponentially distributed)
    fig <- ggplot(bw_distr, aes(x=Betweenness, y=Freq)) + geom_point(colour="grey60")
    fig <- fig + scale_x_log10() + scale_y_log10()
    return (fig)
}


hk_betweenness_test <- function(ppi_name = "string", expr_name="gene_atlas", threshold=0.1)
{
    # get gene properties together with their expression and total counts
    data <- get_node_properties_exprcount(ppi_name, expr_name)

    # categorize into hk and nonHK based on the threshold
    nTissues <- max(data$TotalCount)
    hk_tissues <- (1-threshold)*nTissues
    # bin into HK and nonHK
    hk_data <- data[which(data$ExpressedCount >= hk_tissues),]
    nonhk_data <- data[which(data$ExpressedCount < hk_tissues),]

    # This analysis assumes a normal distribution, which does not fit the
    # actual data well. We still perform this analysis to show the results of
    # Lin et al. (2009)

    # global model:
    my <- mean(nonhk_data$Betweenness)
    s  <- sd(nonhk_data$Betweenness)

    # get z_scores
    z_scores <- (hk_data$Betweenness - my) / s

    # mean z_score (TODO: is this what the Lin paper does?? this would be soo
    # wrong)
    avg_z_score <- mean(z_scores)

    # FIXME:
    # this is a dumb test, a two sample t-test is more appropriate (assuming
    # normal distribution), but data isn't even normally distributed

    return(avg_z_score)
}


hk_ts_betweenness_test <- function(ppi_name = "string", expr_name="gene_atlas", threshold=0.1)
{
    # get gene properties together with their expression and total counts
    data <- get_node_properties_exprcount(ppi_name, expr_name)

    # categorize into hk and nonHK based on the threshold
    nTissues <- max(data$TotalCount)
    hk_tissues <- (1-threshold)*nTissues

    all_tissues <- unique(data$Tissue)

    for (t in all_tissues)
    {
        # filter only the tissue we will look at
        ts_data <- data[which(data$Tissue == t),]
        # bin into HK and nonHK
        hk_data <- ts_data[which(ts_data$ExpressedCount >= hk_tissues),]
        nonhk_data <- ts_data[which(ts_data$ExpressedCount < hk_tissues),]
        # global model:
        my <- mean(nonhk_data$Betweenness)
        s  <- sd(nonhk_data$Betweenness)

        # get z_scores
        z_scores <- (hk_data$Betweenness - my) / s

        # mean z_score (TODO: is this what the Lin paper does?? this would be soo
        # wrong)
        avg_z_score <- mean(z_scores)

        # TODO pass out the scores
    }
}



