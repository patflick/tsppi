
# load required libs:
library(ggplot2)
library(plyr) # for `round_any`
library(gridExtra) # for `grid.arrange`
library(Hmisc) # for spearman correlation + p-value
library(reshape2) # for `melt`

# returns a list of all PPIs that can be analyzed
# TODO: put this as a unified function somewhere (using the SQL database)
# TODO: also so that the order of PPIs is always the same for all plots
get_ppis <- function()
{
    # FIXME: do this properly (not hardcoded)
    ppis <- c("bossi", "string", "psicquic_all", "havu", "ccsb")
    return(ppis)
}

get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all")
    return(exprs)
}

###########################################
#  functions for calculating error bars:  #
###########################################


# standard error:
std <- function(x) sd(x)/sqrt(length(x))

# mean plus/minus one standard error
upper_se <- function(x) mean(x) + std(x)
lower_se <- function(x) mean(x) - std(x)



##########################
#  Bossi & Lehner Fig1B  #
##########################

# data getter:
get_expression_property_data <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="degree_test")
{

    # load sql config and get connection
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # create table names
    prop_table <- paste("prop", ppi_name, expr_name, prop_name, sep="_")
    expr_count_table <- paste(expr_name,"core_expr_counts", sep="_")

    # create SQL query to get ExpressionCount vs Degree
    query <- paste("SELECT a.Gene, a.Value, b.ExpressedCount, b.TotalCount
                    FROM ", prop_table, " AS a",
                   "INNER JOIN ", expr_count_table, " AS b",
                   "ON a.Gene = b.Gene WHERE ExpressedCount > 0")

    # load ppi network from db
    data <- dbGetQuery(con, query)
    return(data);
}


# reproducing the Bossi&Lehner plot (Fig. 1 B)
plot_degree_vs_expr <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="maxts_degree")
{

    data <- get_expression_property_data(ppi_name, expr_name, prop_name)

    # round off the exressedcount as simple "binning"
    # TODO: might be able to do this better :)
    max_expr_count <- max(data$ExpressedCount)
    round_to <- 10
    if (max_expr_count < 60)
    {
        round_to <- round(max_expr_count / 7)
    }
    # round off on x-axis to create bins
    # TODO: this might be possible with ggplot
    data$ExpressedCount <- round_any(data$ExpressedCount, round_to)

    error_width <- max_expr_count / 20

    print(paste(ppi_name, expr_name,":",length(data$ExpressedCount)))

    # plot using ggplot2
    fig <- ggplot(data, aes(x=ExpressedCount, y=Value)) +
          stat_summary(fun.y="mean", geom="line") +
          stat_summary(fun.y="mean", geom="point") +
          # draw error bars:
          stat_summary(fun.ymin="lower_se", fun.ymax="upper_se",
                       geom="errorbar", width=error_width) +
          labs(title=paste("PPI:", ppi_name, ", Expr:", expr_name)) +
          xlab("Num. Tissues") +
          ylab("Degree")

    return(fig)
}

plot_tissue_expr_count_hist <- function(expr_name="gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    expr_count_table <- paste(expr_name, "core_expr_counts", sep="_")
    #expr_count_table <- paste(expr_name, "expr_counts", sep="_")
    query <- paste("SELECT * FROM ", expr_count_table, " ORDER BY ExpressedCount, TotalCount")

    data <- dbGetQuery(con, query)
    data$Gene <- 1:length(data$Gene)
    data <- melt(data, id.vars="Gene")

    fig <- ggplot(data, aes(x=Gene, y=value, group=variable)) +
            labs(title=expr_name) +
            xlab("Genes (sorted by expression)") +
            ylab("Tissue expression") +
            # plot as line + area underneath
            geom_area(aes(fill=variable),position="identity") +
            scale_fill_manual(values=c("gray70", "gray20"))
    return(fig)
}


# Main reason for this plot: to show why HPA/HPA_All shows so few interactions
# between TS and HK
plot_min_max_by_threshold <- function(ppi_name="bossi", expr_name="gene_atlas")
{
    data <- get_min_max_neighbor_expr_data(ppi_name, expr_name)

    nTissues <- max(data$TotalCount)
    threshold <- 0.20
    nTS <- sum(data$ExpressedCount <= threshold * nTissues & data$ExpressedCount != 0)
    n_ts_hk_neighbor <- sum(data$ExpressedCount <= threshold * nTissues& data$ExpressedCount != 0 & data$Max > (1-threshold)*nTissues)

    fig <- ggplot(data, aes(x=ExpressedCount, group=ExpressedCount))# + geom_ribbon(aes(ymin=min(data$Max),ymax=max(data$Max)))
    return (fig)
}

# checking for bossi results:
# 1.) tissue specific proteins make fewer protein interactions than widely
#     expressed proteins (use spearman)
# 2.) most recently evolved proteins have fewer interactions than ancient proteins TODO
test_bossi_1 <- function(ppi_name="bossi", expr_name="gene_atlas", prop_name="coexpr_degree")
{
    data <- get_expression_property_data(ppi_name, expr_name, prop_name)

    # run a spearman test:

    #for (test in c("spearman", "kendall", "pearson"))
    #{
    #    print(cor.test(data$ExpressedCount, data$Value, method=test))
        # TODO plot/output results in publishable format
        #rho = tst$estimate
        #pvalue = tst$p.value
        #print(paste("runing", test, "test: ", ppi_name, expr_name, "rho=", rho, ", p-value=", pvalue))
    #}
    tst <- cor.test(data$ExpressedCount, data$Value, method="spearman")
    #return(tst)
    return(list(value=tst$p.value, label=paste("rho =", round(tst$estimate,2), "\np ~",sprintf("%.1e",tst$p.value))))
}

test_bossi_1_ttest <- function(ppi_name="bossi", expr_name="gene_atlas", prop_name="coexpr_degree")
{
    data <- get_expression_property_data(ppi_name, expr_name, prop_name)


    #max_expr_count <- max(data$ExpressedCount)
    nTissues <- max(data$TotalCount)
    threshold <- 0.2
    # bin into TS and HK
    ts <- data$Value[which(data$ExpressedCount <= threshold*nTissues)]
    hk <- data$Value[which(data$ExpressedCount >= (1-threshold)*nTissues)]

    tst <- t.test(ts,hk)
    #return(tst)
    return(list(value=tst$p.value, label=paste("p ~",sprintf("%.1e",tst$p.value))))
}


get_hkts_sym_threshold_table <- function(expr_name="gene_atlas", threshold=0.1)
{
    expr_count_table <- paste(expr_name, "expr_counts", sep='_')
    # return top and bottom threshold percent as HK and TS (respectively)
    hkts_table <- paste(
                   "SELECT Gene, CASE ",
                   "WHEN ExpressedCount*1.0/TotalCount <= ", threshold,
                   " THEN 'TS'",
                   "WHEN ExpressedCount*1.0/TotalCount >= 1.0 - ", threshold,
                   " THEN 'HK'",
                   "ELSE 'None' END AS Class",
                   "FROM ", expr_count_table)
    hkts_table <- paste("(", hkts_table, ")")
    return (hkts_table)
}


get_min_max_neighbor_expr_data <- function(ppi_name="bossi", expr_name="gene_atlas")
{

    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # have to load from `core`, as the min and max counts are also
    # calculated in the `core` of the expression data set.
    expr_count_table <- paste(expr_name, "core_expr_counts", sep="_")
    min_table <- paste("prop", ppi_name, expr_name, "min_neighbor_expr_count", sep="_")
    max_table <- paste("prop", ppi_name, expr_name, "max_neighbor_expr_count", sep="_")
    # Joining ExpressionCount with Min and Max
    # Min and Max are the minimum and maximum tissue expression of all neighbors
    # (i.e. the minimum and maximum amount of tissues a neighbor of the
    # according gene is expressed in)
    query <- paste("SELECT a.Gene, a.ExpressedCount, a.TotalCount, ",
                   "b.Value AS Min, c.Value AS Max",
                   "FROM ", expr_count_table, " AS a INNER JOIN ",
                   min_table, " AS b ON a.Gene = b.Gene INNER JOIN ",
                   max_table, " AS c ON a.Gene = c.Gene ",
                   # TODO: this excludes all genes that are not expressed at all
                   #       the TODO here is to figure out if this is the right
                   #       position to do this
                   "WHERE ExpressedCount > 0")
    data <- dbGetQuery(con, query)
    return(data)
}


plot_tiles_for_ppi_expr <- function(plot_data, font_size=4)
{

    # create labels if they don't exist yet
    if (length(plot_data$label) == 0)
    {

        if (max(plot_data$value) <= 1.0)
        {
            plot_data$value <- plot_data$value * 100.0
        }

        # generate labels (adding the sizes to the empty rows/cols
        plot_data$label = paste(round(plot_data$value, 1), " %", sep="")
    }

    # plot as heatmap table
    # TODO generalize the heatmap table into a separate function
    p <- ggplot(plot_data, aes(expr, ppi)) +
            geom_tile(aes(fill=value,label=label,colour=value)) +
            # fill colors:
            scale_fill_gradient2("Percent", midpoint=50, limits=c(0,100)) +
            # font colors:
            scale_color_gradient2(low="grey40", mid="black", high="black", limits=c(0,100), guide='none') +
            # set labels:
            # TODO set size here with parameter
            geom_text(aes(label=label, colour=value), size=font_size) +
            xlab('') + ylab('') +
            theme(# disable grid and axis ticks
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  axis.ticks=element_blank(),
                  panel.background=element_blank(),
                  # flip x axis
                  #axis.text.x=element_text(angle=90, hjust=1)
                  axis.text.x=element_text(size=10)
                  )
    return(p)
}

fill_ppi_expr_dataframe <- function(func, ...)
{
    ppis <- get_ppis()
    exprs <- get_exprs()
    len <- length(ppis)*length(exprs)
    data_frame <- data.frame(ppi=rep("",len), expr=rep("",len),value=as.double(rep(0.0,len)),stringsAsFactors=FALSE)

    i <- 1
    for (p in ppis)
    {
        for (e in exprs)
        {
            data_frame$ppi[i] <- p
            data_frame$expr[i] <- e

            result <- func(p,e, ...)
            if (length(result) == 1)
            {
                # assume this is only the value (no label)
                data_frame$value[i] = as.double(result)
            } else {
                # get named variables from list
                data_frame$value[i] <- result$value
                data_frame$label[i] <- result$label
            }

            i <- i + 1
        }
    }
    return(data_frame)
}


# finds the TS/HK classification threshold for which at least 50% of tissue
# specific proteins interact with housekeeping proteins
where_reaches_50_percent_interaction <- function(ppi_name="bossi", expr_name="gene_atlas")
{
    # get the minimum and maximum neighbor tissue expression values
    data <- get_min_max_neighbor_expr_data(ppi_name, expr_name)

    nTissues <- max(data$TotalCount)

    target <- 0.5

    # find the appropriate threshold by binary search
    # FIXME: don't need binary search, there's only a fixed (and small)
    #        number of tissues. I can just iterate through them linearly
    left <- 0.05
    right <- 0.5

    while (right - left > 0.01)
    {
        threshold <- (right+left)/2.0
        nTS <- sum(data$ExpressedCount <= threshold * nTissues)
        n_ts_hk_neighbor <- sum(data$ExpressedCount <= threshold * nTissues & data$Max > (1-threshold)*nTissues)

        value <- n_ts_hk_neighbor / nTS

        if (value > target)
        {
            # go in lower half
            right <- threshold
        } else {
            left <- threshold
        }
    }
    print(paste(ppi_name, expr_name, " threshold for 50% TS->HK: ", threshold, " with actual value: ", value))
}


# 2.) to which extend to tissue-specific proteins interact with most widely expressed proteins?

# TODO: ok this can be plotted now, the plot has to be made `nicer`
#       especially considering the graident colors. Also add more than one label
#       per square
test_bossi_2 <- function(ppi_name="bossi", expr_name="gene_atlas", threshold=0.125)
{
    # get the data
    data <- get_min_max_neighbor_expr_data(ppi_name, expr_name)

    nTissues <- max(data$TotalCount)

    nTS <- sum(data$ExpressedCount <= threshold * nTissues & data$ExpressedCount != 0)
    n_ts_hk_neighbor <- sum(data$ExpressedCount <= threshold * nTissues &
                            data$ExpressedCount != 0 &
                            data$Max > (1-threshold) * nTissues)

    # TODO actually output/plot this data
    print(paste(ppi_name, expr_name, ": ", n_ts_hk_neighbor, "/", nTS, " max expr neighbor ", max(data$Max), "/",nTissues, " total size:", length(data$ExpressedCount) ))
    result <- n_ts_hk_neighbor*100.0/nTS
    return (list(value=result, label=paste(round(result, 1), "%")))
}


test_bossi_2_expectation <- function(ppi_name="bossi", expr_name="gene_atlas", threshold=0.125, account_for_ts=FALSE)
{
    # get the data
    # TODO: test different degree sources
    if (account_for_ts){
        data <- get_expression_property_data(ppi_name, expr_name, "coexpr_degree")
    } else {
        data <- get_expression_property_data(ppi_name, expr_name, "degree")
    }
    data <- data[which(data$ExpressedCount != 0),]

    nTissues <- max(data$TotalCount)

    nTS <- sum(data$ExpressedCount <= threshold * nTissues)
    nHK <- sum(data$ExpressedCount >= (1-threshold) * nTissues)
    N <- length(data$ExpressedCount)
    mean_degree <- mean(data$Value)

    # probability of an edge to HK
    # ----------------------------
    # use the hypergeometric distribution to calc the expectation for an
    # edge of a TS node to connect to HK
    # thus we want the probability P(X >= 1) using the hyper geometric distr
    # and:
    # P(X >= 1) = 1 - P(X = 0)
    data$HkEdgeProb <- (1 - dhyper(0, nHK, (N-1) - nHK, data$Value))

    if (account_for_ts) {
        # take TS smaller degrees into account
        ex_edge_to_hk_exists <- mean(data$HkEdgeProb[which(data$ExpressedCount <= threshold * nTissues)])
    } else {
        # average over all expectations
        ex_edge_to_hk_exists <- mean(data$HkEdgeProb)
    }

    # now multiply by the number of tissue specific nodes, in order to
    # get the expected number of TS having an edge to HK
    ex_n_ts_hk_neighbor <- ex_edge_to_hk_exists * nTS

    # WRONG!
    # ex_n_ts_hk_neighbor <- mean_degree * (nHK*1.0/N) * nTS

    # TODO: test this function
    result <- ex_n_ts_hk_neighbor*100.0/nTS
    return (list(value=result, label=paste(round(result, 1), "%")))
}


# testing bossi 3:
#   almost all houskeeping genes have some non housekeeping interaction
test_bossi_3 <- function(ppi_name="bossi", expr_name="gene_atlas", threshold=0.125)
{
    # get the data
    data <- get_min_max_neighbor_expr_data(ppi_name, expr_name)

    nTissues <- max(data$TotalCount)

    nHK <- sum(data$ExpressedCount >= (1-threshold) * nTissues & data$ExpressedCount != 0)
    # get the number of HK proteins that interact directly with non HK proteins
    n_hk_non_hk_neighbor <- sum(data$ExpressedCount >= (1-threshold) * nTissues &
                                data$ExpressedCount != 0 &
                                data$Min < (1-threshold) * nTissues)

    # get result (percentage of HK that interact with non HK)
    result <- n_hk_non_hk_neighbor*100.0/nHK
    return (list(value=result, label=paste(round(result, 1), "%")))
}

test_bossi_3_expectation <- function(ppi_name="bossi", expr_name="gene_atlas", threshold=0.125, account_for_hk=FALSE)
{
    # get the data
    if (account_for_hk){
        data <- get_expression_property_data(ppi_name, expr_name, "coexpr_degree")
    } else {
        data <- get_expression_property_data(ppi_name, expr_name, "degree")
    }
    data <- data[which(data$ExpressedCount != 0),]

    nTissues <- max(data$TotalCount)

    nTS <- sum(data$ExpressedCount <= threshold * nTissues)
    nHK <- sum(data$ExpressedCount >= (1-threshold) * nTissues)
    N <- length(data$ExpressedCount)
    mean_degree <- mean(data$Value)

    # probability of an edge to non-HK
    # ----------------------------
    # use the hypergeometric distribution to calc the expectation for an
    # edge of a HK node to connect to non-HK
    # thus we want the probability P(X >= 1) using the hyper geometric distr
    # and:
    # P(X >= 1) = 1 - P(X = 0)
    data$NonHkEdgeProb <- 1 - dhyper(0, N - nHK, nHK, data$Value)

    if (account_for_hk) {
        # take the degree distribution of HK proteins into account
        ex_edge_to_nonhk_exists <- mean(data$NonHkEdgeProb[which(data$ExpressedCount >= (1 - threshold) * nTissues)])
    } else {
        # average over all expectations
        ex_edge_to_nonhk_exists <- mean(data$NonHkEdgeProb)
    }

    # now multiply by the number of tissue specific nodes, in order to
    # get the expected number of TS having an edge to HK
    ex_n_hk_nonhk_neighbor <- ex_edge_to_nonhk_exists * nHK

    # TODO: test this function
    result <- ex_n_hk_nonhk_neighbor*100.0/nHK
    return (list(value=result, label=paste(round(result, 1), "%")))
}




get_hkts_edge_summary <- function(ppi_name="bossi", expr_name="gene_atlas", prop_name="maxts_degree")
{
    # this can be solved for various definitions of "housekeeping" and
    # tissue specific proteins by only using SQL:

    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # what's the definition of tissuespecific vs. housekeeping?
    hkts_table <- get_hkts_sym_threshold_table(expr_name, 0.2)


    # TODO (test): min(Class1,Class2),max(Class1,Class2)
    query <- paste("SELECT MIN(a.Class,b.Class) AS Class1, MAX(a.Class, b.Class) AS Class2,",
                   "COUNT(*) AS Count FROM ",
                   ppi_name, " AS c INNER JOIN ",
                   hkts_table, " AS a ON c.Gene1 = a.Gene INNER JOIN ",
                   hkts_table, " AS b ON c.Gene2 = b.Gene",
                   "GROUP BY Class1, Class2")
    data <- dbGetQuery(con, query);
    return (data)
}


create_test_degrees <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="degree_test")
{

    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # create table names
    prop_table <- paste("prop", ppi_name, expr_name, prop_name, sep="_")

    dbSendQuery(con, paste("DROP TABLE IF EXISTS ", prop_table))
    query <- paste("CREATE TABLE IF NOT EXISTS ", prop_table, " AS ",
                   "SELECT Gene, COUNT(*) AS Value FROM ",
                   "(",
                   "SELECT Gene1 AS Gene FROM ", ppi_name,
                   # requires UNION all, otherwise all degrees will be 1
                   " UNION ALL ",
                   "SELECT Gene2 As Gene FROM ", ppi_name,
                   ")",
                   "GROUP BY Gene")
    dbSendQuery(con, query);
}


compare_test_degrees <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="degree_test")
{

    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # create table names
    prop_table_1 <- paste("prop", ppi_name, expr_name, "_degree", sep="_")
    prop_table_2 <- paste("prop", ppi_name, expr_name, "degree_test", sep="_")

    query <- paste("SELECT a.Gene, a.Value AS Value1, b.Value AS Value2 ",
                   "FROM ", prop_table_1, " AS a INNER JOIN ",
                   prop_table_2, " AS b ON a.Gene = b.Gene")
    data <- dbGetQuery(con, query)

    #print(which(data$Value1 != data$Value2))
}

test_for_all <- function(func)
{
    ppis <- get_ppis()
    exprs <- get_exprs()
    for (p in ppis)
    {
        for (e in exprs)
        {
            func(p,e)
        }
    }
}

for_all_expr <- function(func)
{
    exprs <- get_exprs()
    for (e in exprs)
    {
        func(e)
    }
}

plot_all_expr <- function(plot_func)
{
    exprs <- get_exprs()
    figs <- list()
    for (e in exprs)
    {
        fig <- plot_func(e)
        figs <- c(figs, list(fig))
    }
    all_figs <- do.call(grid.arrange, figs)
    return(all_figs)
}

plot_all <- function(plot_func,...)
{
    ppis <- get_ppis()
    exprs <- get_exprs()

    par(mfrow=c(length(ppis),length(exprs)))

    figs <- list()
    for (p in ppis)
    {
        for (e in exprs)
        {
            # call the bossi plotting function
            #fig <- plot_degree_vs_expr(p, e, "coexpr_degree")
            fig <- plot_func(p,e,...)
            #fig <- plot_degree_vs_expr(p, e, "degree_test")
            figs <- c(figs, list(fig))
        }
    }

    #all_figs <- do.call(grid.arrange, figs)
    all_figs <- do.call(arrangeGrob, figs)

    return(all_figs)
    # in order to pass named parameters to grid.arrange:
    #do.call(grid.arrange, c(figs,list(nrow=5,left="hello there"))
}

#######################################################################
#                         random expectation                          #
#######################################################################


plot_expectation_corr <- function(data, threshold)
{
    # get linear regression model
    l <- lm(expected_value~value, data=data)
    lm_intercept <- l$coefficients[1]
    lm_slope <- l$coefficients[2]

    r <- cor(data$value, data$expected_value)

    fig <- ggplot(data, aes(x = value, y = expected_value)) +
        geom_point() +
        coord_fixed(1, xlim=c(0,100), ylim=c(0,100)) +
        geom_abline(intercept = 0, slope = 1, colour="grey60",lty="dotted") +
        geom_abline(intercept = lm_intercept, slope = lm_slope, colour="grey40") +
        annotate("text", x = 20, y = lm_intercept + 30*lm_slope, label=paste("r =",round(r,3)), size=3.5) +
        xlab("Actual") +
        ylab("Random expectation") +
        labs(title=paste("Expected vs real interactions (t =", 100*threshold, "%)"))

    return(fig)
}

bossi_2_expectation_corr <- function(threshold = 0.125, account_for_ts=FALSE)
{
    data_o <- fill_ppi_expr_dataframe(test_bossi_2, threshold)
    data_e <- fill_ppi_expr_dataframe(test_bossi_2_expectation, threshold, account_for_ts)

    data <- data_o
    data$expected_value <- data_e$value

    return(data)
}

bossi_3_expectation_corr <- function(threshold = 0.125, account_for_hk=FALSE)
{
    data_o <- fill_ppi_expr_dataframe(test_bossi_3, threshold)
    data_e <- fill_ppi_expr_dataframe(test_bossi_3_expectation, threshold, account_for_hk)

    data <- data_o
    data$expected_value <- data_e$value

    return(data)
}

#######################################################################
#                 Generate actual plots (into PDF)                    #
#######################################################################

bossi_1 <- function()
{
    # TODO: put all properties into one plot..
    for (prop in c("maxts_degree", "coexpr_degree","degree"))
    {
        p <- plot_all(plot_degree_vs_expr, prop)
        # TODO size approriately
        pdf(paste("../figs/bossi1_all_prop_",prop,".pdf",sep=""), width=20, height=15)
        print(p)
        dev.off()
    }

    # output small examples
    # TODO: what would the random expectation look like?
    p <- plot_degree_vs_expr("bossi", "gene_atlas", "coexpr_degree")
    pdf("../figs/bossi1_bossi_gene_atlas.pdf", width=5, height=3)
    p <- p + labs(title="Interaction degree by tissue expression")
    p <- p + xlab("Number of tissues in which protein is expressed")
    p <- p + ylab("Mean interaction degree")
    print(p)
    dev.off()
    # again with other data
    p <- plot_degree_vs_expr("ccsb", "hpa", "coexpr_degree")
    pdf("../figs/bossi1_ccsb_hpa.pdf", width=5, height=3)
    p <- p + labs(title="Interaction degree by tissue expression")
    p <- p + xlab("Number of tissues in which protein is expressed")
    p <- p + ylab("Mean interaction degree")
    print(p)
    dev.off()
    # for string <-> emtab
    p <- plot_degree_vs_expr("string", "emtab", "coexpr_degree")
    pdf("../figs/bossi1_string_emtab.pdf", width=5, height=3)
    p <- p + labs(title="Interaction degree by tissue expression")
    p <- p + xlab("Number of tissues in which protein is expressed")
    p <- p + ylab("Mean interaction degree")
    print(p)
    dev.off()


    # spearmans rho for all PxE and for different degree calculations
    for (prop in c("maxts_degree", "coexpr_degree","degree"))
    {
        pdf(paste("../figs/bossi1_rho_", prop, ".pdf", sep=""), width=7, height=3.5)
        data <- fill_ppi_expr_dataframe(test_bossi_1,prop)
        p <- plot_tiles_for_ppi_expr(data, 3)
        p <- p + scale_fill_gradient2("P-value",na.value="white",midpoint=1, mid="red",high="red", low="white",limits=c(1e-3,1e-1),trans="log10")
        p <- p + labs(title="Correlation test: Spearman's rho")
        plot(p)
        dev.off()
    }

    # simple t-test for TS vs. HK degrees
    for (prop in c("maxts_degree", "coexpr_degree","degree"))
    {
        pdf(paste("../figs/bossi1_ttest_", prop, ".pdf", sep=""), width=7, height=3)
        data <- fill_ppi_expr_dataframe(test_bossi_1_ttest, prop)
        p <- plot_tiles_for_ppi_expr(data, 3)
        p <- p + scale_fill_gradient2("P-value",na.value="white",midpoint=1, mid="red",high="red", low="white",limits=c(1e-3,1e-1),trans="log10")
        p <- p + labs(title="T-test: degrees of TS vs. HK proteins")
        plot(p)
        dev.off()
    }
}


bossi_2 <- function()
{
    for (t in c(0.1, 0.125, 0.15, 0.2, 0.5))
    {
        data <- fill_ppi_expr_dataframe(test_bossi_2, t)
        p <- plot_tiles_for_ppi_expr(data)
        p <- p + labs(title=paste("TS that interact with HK (threshold:", t*100,"%)"))

        pdf(paste("../figs/bossi2_all_t",t*100,".pdf",sep=""), width=7, height=3)
        print(p)
        dev.off()
    }
    # TODO? output single plot examples

    # Random expectation:
    for (t in c(0.1, 0.125, 0.15, 0.2, 0.5))
    {
        # random expectation percentage tiles
        data <- fill_ppi_expr_dataframe(test_bossi_2_expectation, t)
        p <- plot_tiles_for_ppi_expr(data)
        p <- p + labs(title=paste("Expected TS that interact with HK (threshold:", t*100,"%)"))

        pdf(paste("../figs/bossi2_all_expected_t",t*100,".pdf",sep=""), width=7, height=3)
        print(p)
        dev.off()

        # random expectation correlation
        data <- bossi_2_expectation_corr(t, FALSE)
        p <- plot_expectation_corr(data, t)
        pdf(paste("../figs/bossi2_expected_corr_t",t*100,".pdf",sep=""), width=4, height=4)
        print(p)
        dev.off()

        # random expectation correlation accounting for TS degree distr
        data <- bossi_2_expectation_corr(t, TRUE)
        p <- plot_expectation_corr(data, t)
        pdf(paste("../figs/bossi2_expected_corr_tsdegree_t",t*100,".pdf",sep=""), width=4, height=4)
        print(p)
        dev.off()
    }
}


bossi_3 <- function()
{
    for (t in c(0.1, 0.125, 0.15, 0.2, 0.5))
    {
        data <- fill_ppi_expr_dataframe(test_bossi_3, t)
        p <- plot_tiles_for_ppi_expr(data)
        p <- p + labs(title=paste("HK that interact with non HK (threshold:", t*100,"%)"))

        pdf(paste("../figs/bossi3_all_t",t*100,".pdf",sep=""), width=7, height=3)
        print(p)
        dev.off()
    }

    # Random expectation
    for (t in c(0.1, 0.125, 0.15, 0.2, 0.5))
    {
        data <- fill_ppi_expr_dataframe(test_bossi_3_expectation, t)
        p <- plot_tiles_for_ppi_expr(data)
        p <- p + labs(title=paste("Expected HK that interact with non HK (threshold:", t*100,"%)"))

        pdf(paste("../figs/bossi3_all_expected_t",t*100,".pdf",sep=""), width=7, height=3)
        print(p)
        dev.off()

        # random expectation correlation
        data <- bossi_3_expectation_corr(t, FALSE)
        p <- plot_expectation_corr(data, t)
        pdf(paste("../figs/bossi3_expected_corr_t",t*100,".pdf",sep=""), width=4, height=4)
        print(p)
        dev.off()

        # random expectation correlation accounting for TS degree distr
        data <- bossi_3_expectation_corr(t, TRUE)
        p <- plot_expectation_corr(data, t)
        pdf(paste("../figs/bossi3_expected_corr_tsdegree_t",t*100,".pdf",sep=""), width=4, height=4)
        print(p)
        dev.off()
    }
}

