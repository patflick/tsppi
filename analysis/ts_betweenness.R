
library(ggplot2)
library(xtable) # for converting a data frame into a latex table
library(gridExtra) # for `grid.arrange`
library(scales) # for muted


source("ppi_utils.R")
source("expr_utils.R")
source("plot_tiles.R")

get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all")
    return(exprs)
}


# returns node properties for ts graphs
get_node_properties <- function(ppi_name = "string", expr_name = "gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn()

    ppi_node_prop_tbl <- paste(ppi_name,expr_name, "node_properties", sep="_")
    query <- paste("SELECT * FROM ", ppi_node_prop_tbl)
    node_props <- dbGetQuery(con, query)

    return(node_props)
}


# get tissue specific node properties
get_ts_node_properties <- function(ppi_name = "string", expr_name = "gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn()

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
    con <- get_sql_conn()

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
    con <- get_sql_conn()

    # create table names
    ppi_node_prop_tbl <- paste(ppi_name,expr_name, "ts_node_properties", sep="_")
    # for TS<->HK classification
    expr_count_table <- paste(expr_name,"core_expr_counts", sep="_")
    # for filtering of expressed genes
    expr_core_table <- paste(expr_name, "core", sep="_")

    # create SQL query to
    #         get only those that are actually expressed in the tissue we are
    #         looking at, otherwise the results will be highly skewedby
    #         including all the degree 0 nodes with betweenness zero
    query <- paste("SELECT a.*, b.ExpressedCount, b.TotalCount FROM ",
                   ppi_node_prop_tbl, " AS a INNER JOIN ", expr_count_table,
                   " AS b ON a.Gene = b.Gene INNER JOIN ", expr_core_table,
                   " AS c ON a.Gene = c.Gene AND a.Tissue = c.Type ",
                   " WHERE c.Expressed = 1")

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


# the test for z-score that Lin et al. use for significance testing of
# HK gene properties within the whole population of genes
# params:
# - x:              a vector of properties (i.e. betweenness)
# - subset_which:   indeces into `x` which define the subset to be tested for
#                   significant differences
# - n_samples:      The number of samplings from the population
subset_significance_lin_et_al <- function(x, subset_which, n_samples=1000)
{
    # This analysis assumes a normal distribution, which does not fit the
    # actual data well. We still perform this analysis to show the results of
    # Lin et al. (2009)
    # Note: the following are the steps performed by Lin et al.:

    # 1) get mean of HK (or subset)
    subset_mean <- mean(x[subset_which])

    # 2a) get number of input values
    n <- length(x)
    # 2b) get the number of HK proteins (-> the subset)
    subset_size <- length(subset_which)

    # 3) sample the whole data 1000 times to get distribution of means
    sampled_means <- rep(0.0, n_samples)
    for (i in 1:n_samples)
    {
        sampled_data <- x[sample(1:n, subset_size)]
        sample_mean <- mean(sampled_data)
        sampled_means[i] <- sample_mean
    }

    # 4) get model parameters of sample means
    my <- mean(sampled_means)
    sigma <- sd(sampled_means)

    # 5) get z-score of HK mean in this model
    z_score <- (subset_mean - my) / sigma

    return (z_score)
}


hk_betweenness_test <- function(ppi_name = "string", expr_name="gene_atlas", threshold=0.1, HK=TRUE)
{
    # get gene properties together with their expression and total counts
    data <- get_node_properties_exprcount(ppi_name, expr_name)

    # categorize into hk and nonHK based on the threshold
    nTissues <- max(data$TotalCount)
    hk_tissues <- (1-threshold)*nTissues
    if (HK)
    {
        which_proteins <- which(data$ExpressedCount >= hk_tissues)
    } else {
        which_proteins <- which(data$ExpressedCount <= threshold*nTissues)
    }

    # call Lin et al. testing
    z_score <- subset_significance_lin_et_al(data$Betweenness, which_proteins)

    # global model:
    #my <- mean(nonhk_data$Betweenness)
    #s  <- sd(nonhk_data$Betweenness)

    # get z_scores
    #z_scores <- (hk_data$Betweenness - my) / s

    # mean z_score (TODO: is this what the Lin paper does?? this would be soo
    # wrong)
    #avg_z_score <- mean(z_scores)

    # FIXME:
    # this is a dumb test, a two sample t-test is more appropriate (assuming
    # normal distribution), but data isn't even normally distributed

    return(z_score)
}


hk_ts_betweenness_test <- function(ppi_name = "string", expr_name="gene_atlas", threshold=0.1, HK=TRUE)
{
    # get gene properties together with their expression and total counts
    data <- get_ts_node_properties_exprcount(ppi_name, expr_name)

    # categorize into hk and nonHK based on the threshold
    nTissues <- max(data$TotalCount)
    hk_tissues <- (1-threshold)*nTissues

    all_tissues <- unique(data$Tissue)
    z_scores <- rep(0.0, length(all_tissues))
    for (i in 1:length(all_tissues))
    {
        t <- all_tissues[i]
        # filter only the tissue we will look at
        ts_data <- data[which(data$Tissue == t),]
        # bin into HK and nonHK
        if (HK){
            hk_which <- which(ts_data$ExpressedCount >= hk_tissues)
        } else {
            hk_which <- which(ts_data$ExpressedCount <= threshold*nTissues)
        }

        # get z_score
        z_score <- subset_significance_lin_et_al(ts_data$Betweenness, hk_which)
        z_scores[i] <- z_score
    }

    return (z_scores)
}


hk_all_scores <- function(threshold=0.1, HK=TRUE)
{
    # prepare the output data frame
    df <- data.frame(ppi=character(0), expr=character(0), value=double(0), stringsAsFactors=FALSE)
    for (p in get_ppis())
    {
        for (e in get_exprs())
        {
            # get the z-scores
            z_score <- hk_betweenness_test(p, e, threshold, HK)
            # add row to data frame
            new_row <- data.frame(ppi=p, expr=e, value=z_score,stringsAsFactors=FALSE)
            df <- rbind(df, new_row)
        }
    }
    return (df)
}

hk_ts_all_scores <- function(threshold=0.1, HK=TRUE, agg_func=NA)
{
    # prepare the output data frame
    df <- data.frame(ppi=character(0), expr=character(0), value=double(0), stringsAsFactors=FALSE)
    for (p in get_ppis())
    {
        for (e in get_exprs())
        {
            # get the z-scores
            z_scores <- hk_ts_betweenness_test(p, e, threshold, HK)
            n <- length(z_scores)
            # add them into the dataset
            if (is.na(agg_func))
            {
                # no aggregation
                new_row <- data.frame(ppi=rep(p,n), expr=rep(e,n), value=z_scores,stringsAsFactors=FALSE)
                df <- rbind(df, new_row)
            } else {
                # aggregate the z_values
                agg_z_score <- agg_func(z_scores)
                new_row <- data.frame(ppi=p, expr=e, value=agg_z_score,stringsAsFactors=FALSE)
                df <- rbind(df, new_row)
            }
        }
    }
    return (df)
}


agg_summary <- function(data)
{
    # get min, max and mean
    sum_data <- aggregate(value~ppi+expr, data=data, FUN=min)
    sum_data$min <- sum_data$value

    xmean <- aggregate(value~ppi+expr, data=data, FUN=mean)
    sum_data$mean <- xmean$value

    xmax <- aggregate(value~ppi+expr, data=data, FUN=max)
    sum_data$max <- xmax$value

    sum_data$value <- sum_data$mean

    return (sum_data)
}


# plot the (min,max,median, quantile) of z-scores
plot_all_zvalues <- function(HK=TRUE)
{
    threshold <- 0.1
    data <- hk_ts_all_scores(threshold, HK)

    # summarize data
    sum_data <- agg_summary(data)
    # create label
    sum_data$label <- paste("[",round(sum_data$min, 1),", ",round(sum_data$max, 1), "]\n", "mean = ", round(sum_data$mean, 1), sep="")

    fig <- plot_tiles_for_ppi_expr(sum_data, 3)
    if (HK)
    {
        fig <- fig + scale_fill_gradient2("z-score",na.value=alpha(muted("blue"),0.9),midpoint=0, mid="white",high=alpha(muted("blue"),0.9), low=muted("red"),limits=c(-5,5))
    }else{
        fig <- fig + scale_fill_gradient2("z-score",na.value=alpha(muted("red"),0.9),midpoint=0, mid="white",high=alpha(muted("blue"),0.9), low=muted("red"),limits=c(-5,5))
    }
    #fig <- fig + scale_fill_gradient2("Percent", midpoint=0, limits=c(-5,5))

    return (fig)
}

# plot the (min,max,median, quantile) of z-scores
plot_zvalues <- function(HK=TRUE)
{
    threshold <- 0.1
    data <- hk_all_scores(threshold, HK)

    # create label
    data$label <- round(data$value, 1)

    fig <- plot_tiles_for_ppi_expr(data)
    if (HK)
    {
        fig <- fig + scale_fill_gradient2("z-score",na.value=alpha(muted("blue"),0.9),midpoint=0, mid="white",high=alpha(muted("blue"),0.9), low=muted("red"),limits=c(-5,5))
    }else{
        fig <- fig + scale_fill_gradient2("z-score",na.value=alpha(muted("red"),0.9),midpoint=0, mid="white",high=alpha(muted("blue"),0.9), low=muted("red"),limits=c(-5,5))
    }
    #fig <- fig + scale_fill_gradient2("Percent", midpoint=0, limits=c(-5,5))

    return (fig)
}


lin_between <- function()
{
    # global ts network HK betweenness
    pdf(paste("../figs/lin_betweenness_hk.pdf", sep=""), width=7, height=3)
    p <- plot_zvalues(TRUE)
    p <- p + labs(title="Betweenness of HK genes: z-scores")
    plot(p)
    dev.off()

    # global ts network TS betweenness
    pdf(paste("../figs/lin_betweenness_ts.pdf", sep=""), width=7, height=3)
    p <- plot_zvalues(FALSE)
    p <- p + labs(title="Betweenness of TS genes: z-scores")
    plot(p)
    dev.off()

    # specific ts network HK betweenness
    pdf(paste("../figs/lin_ts_betweenness_hk.pdf", sep=""), width=7, height=3.5)
    p <- plot_all_zvalues(TRUE)
    p <- p + labs(title="Betweenness of HK genes in all subnetworks: z-scores")
    plot(p)
    dev.off()

    # specific ts network TS betweenness
    pdf(paste("../figs/lin_ts_betweenness_ts.pdf", sep=""), width=7, height=3.5)
    p <- plot_all_zvalues(FALSE)
    p <- p + labs(title="Betweenness of TS genes in all subnetworks: z-scores")
    plot(p)
    dev.off()
}
