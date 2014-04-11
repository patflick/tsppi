library(ggplot2)
library(reshape2) # for `melt`
library(gridExtra) # for `grid.arrange`

# load Expr name mapping
source("expr_utils.R")

get_all_expr_size_stats <- function()
{
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    for (e in get_exprs())
    {
        query <- paste("SELECT COUNT(DISTINCT Gene), COUNT(DISTINCT Type) FROM ", e)
        data <- dbGetQuery(con, query)

        nGenes <- data[1,1]
        nTypes <- data[1,2]

        cat(sprintf("%s    & %i    &     %i  \\\\\n", to_expr_name(e), nGenes, nTypes))
    }
}

# get the expression values table for the given expression data set
get_expr_values <- function(expr_name="gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    expr_normalized_table <- paste(expr_name, "normalized", sep="_")
    query <- paste("SELECT * FROM ", expr_normalized_table)

    data <- dbGetQuery(con, query)

    return (data)
}


plot_hist_normalized_expr <- function(expr_name="gene_atlas")
{
    # get the data
    data <- get_expr_values(expr_name)

    #x <- data$ExpressionValue
    #print(paste(expr_name, " median:", quantile(x)))
    #x <- x[which(x <= 5)]
    #hist(x, xlab=expr_name)
    #plot(ecdf(x), main=expr_name)
    fig <- ggplot(data, aes(x = ExpressionValue))
    fig <- fig + stat_ecdf()
    fig <- fig + coord_cartesian(xlim = c(-0.5, quantile(x,c(0.90))[1]))
    #fig <- fig + scale_x_log10()
    fig <- fig + labs(title=expr_name)
    return(fig)
}


plot_tissue_expr_count_hist <- function(expr_name="gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    expr_count_table <- paste(expr_name, "expr_counts", sep="_")
    query <- paste("SELECT * FROM ", expr_count_table, " ORDER BY ExpressedCount*1.0/TotalCount")

    data <- dbGetQuery(con, query)
    data$Gene <- 1:length(data$Gene)

    fig <- ggplot(data, aes(x=Gene, y=ExpressedCount/TotalCount)) +
            labs(title=to_expr_name(expr_name)) +
            xlab("Genes") +
            ylab("Tissue expression") +
            # plot as line + area underneath
            geom_area(position="identity") 
            #scale_fill_manual(values=c("gray70", "gray20"))
    return(fig)
}


plot_tshk_example <- function(expr_name="hpa", threshold=0.15)
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    expr_count_table <- paste(expr_name, "expr_counts", sep="_")
    query <- paste("SELECT * FROM ", expr_count_table, " WHERE ExpressedCount > 0 ORDER BY ExpressedCount*1.0/TotalCount")

    data <- dbGetQuery(con, query)
    # replace gene names with numbers from {1,...,n}
    data$Gene <- 1:length(data$Gene)

    expr_frac <- data$ExpressedCount/data$TotalCount

    # find the index for thresholds t and 1-t
    ts_index <- binsearch_for(expr_frac, threshold)
    hk_index <- binsearch_for(expr_frac, 1-threshold)

    fig <- ggplot(data, aes(x=Gene, y=ExpressedCount/TotalCount)) +
            labs(title=paste("TS vs HK (t = ",100*threshold, "%)",sep="")) +
            xlab("Genes") +
            ylab("Tissue expression") +
            # plot as line + area underneath
            geom_area(position="identity") 

    # add annotations
    t <- threshold
    # TS
    fig <- fig + annotate("segment", x = 0,  xend = ts_index, y = t, yend = t, colour = "blue", lty="dashed")
    fig <- fig + annotate("text", y = t, x=0, label=paste("t = ",t*100,"%",sep=""), vjust=-1, hjust=0, size=3.8,  colour="blue")
    fig <- fig + annotate("segment", x = ts_index,  xend = ts_index, y = -.1, yend = 1, colour = "blue", lty="dashed")
    fig <- fig + annotate("segment", x = ts_index,  xend = 0, y = 0.5, yend = 0.5, colour = "blue", arrow=arrow(type="closed",length = unit(0.1, "inches")))
    fig <- fig + annotate("text", y = 0.5, x=ts_index/2, label="TS", vjust=-1, size=5,  colour="blue")
    # HK
    fig <- fig + annotate("segment", x = 0,  xend = hk_index, y = 1-t, yend = 1-t, colour = "red", lty="dashed")
    fig <- fig + annotate("text", y = 1-t, x=0, label=paste("t = ",(1-t)*100,"%",sep=""), vjust=-1, hjust=0, size=3.8,  colour="red")
    fig <- fig + annotate("segment", x = hk_index,  xend = hk_index, y = -.1, yend = 1, colour = "red", lty="dashed")
    fig <- fig + annotate("segment", x = hk_index,  xend = hk_index+500, y = 0.5, yend = 0.5, colour = "red", arrow=arrow(type="closed",length = unit(0.1, "inches")))
    fig <- fig + annotate("text", y = 0.5, x=hk_index+250, label="HK", vjust=-1, size=5,  colour="red")


    return(fig)
}

plot_tshk_thresholds <- function(expr_name="hpa")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    expr_count_table <- paste(expr_name, "expr_counts", sep="_")
    query <- paste("SELECT * FROM ", expr_count_table, " WHERE ExpressedCount > 0 ORDER BY ExpressedCount*1.0/TotalCount")

    data <- dbGetQuery(con, query)
    # replace gene names with numbers from {1,...,n}
    data$Gene <- 1:length(data$Gene)

    expr_frac <- data$ExpressedCount/data$TotalCount

    # get number of total genes
    query <- paste("SELECT COUNT(DISTINCT Gene) FROM ", expr_count_table)
    cdata <- dbGetQuery(con, query)
    nGenes <- cdata[1,1]

    # get the the TS threshold -> size data
    ts_x <- unique(expr_frac[which(expr_frac <= 0.5)])
    ts_y <- ts_x
    for (i in 1:length(ts_x))
    {
        ts_y[i] <- binsearch_for(expr_frac, ts_x[i]) / nGenes
    }

    # get the the HK threshold -> size data
    hk_x <- 1.0 - unique(expr_frac[which(expr_frac >= 0.5)])
    hk_y <- hk_x
    max_index = length(data$Gene)
    for (i in 1:length(ts_x))
    {
        hk_y[i] <- (max_index - binsearch_for(expr_frac, 1- hk_x[i])) / nGenes
    }

    # combine TS and HK graphs into one dataframe for ggplot
    df_ts <- data.frame(Class=rep("TS", length(ts_x)), x=ts_x, y=ts_y)
    df_hk <- data.frame(Class=rep("HK", length(hk_x)), x=hk_x, y=hk_y)

    df <- rbind(df_ts, df_hk)


    fig <- ggplot(df, aes(x=x, y=y, group=Class,colour=Class,lty=Class)) +
            geom_line() +
            labs(title=to_expr_name(expr_name)) +
            xlab("Threshold") +
            ylab("Rel. size")
            # plot as line + area underneath

    return(fig)
}

# gets the discrete cds of the data
inverse_cdf <- function(data)
{
    # first sort the data and get the unique values
    data <- sort(data)
    udata <- sort(unique(data))
    freq <- udata
    n <- length(data)


    j <- 1
    for (i in 1:length(udata))
    {
        # find the first position in `data` to be equal to this
        while (data[j] != udata[i])
        {
            j <- j + 1
        }
        # now find the last position to still equal this
        while (j < n && data[j] == udata[i])
        {
            j <- j + 1
        }
        # set the cdf freq
        freq[i] = (j-1) / n
    }

    # return the dataframe for the cdf
    df <- data.frame(val=udata, freq=freq)
    return (df)
}

binsearch_for <- function(data, value)
{
    n <- length(data)
    # binary search
    left <- 1
    right <- n
    while (left < right - 1)
    {
        # integer divide
        mid <- (right + left)%/%2
        if (data[mid] < value)
        {
            # bottom half
            left <- mid
        }
        else
        {
            # top half
            right <- mid
        }
    }

    # return the index of the first element with the given value
    if (data[left] == value)
    {
        return(left)
    }
    else
    {
        return(right)
    }
}

get_thres_freq <- function(df, t)
{
    n <- length(df$val)
    # binary search
    left <- 1
    right <- n
    while (left < right - 1)
    {
        # integer divide
        mid <- (right + left)%/%2
        if (df$val[mid] < t)
        {
            # bottom half
            left <- mid
        }
        else
        {
            # top half
            right <- mid
        }
    }
    return (df$freq[left])
}

# plots the percentage of expressed classifications linearly on the x-axis
plot_expr_inv_cdf <- function(expr_name="gene_atlas", threshold = 100, t_lbl=threshold)
{
    # get the data
    data <- get_expr_values(expr_name)

    #y <- unique(data$ExpressionValue)
    #
    #inv_cdf <- function(val)
    #{
    #    return (sum(data$ExpressionValue <= val) / length(data$ExpressionValue))
    #}

    cdf_data <- inverse_cdf(data$ExpressionValue)

    # get the freq of the cdf at the threshold position
    thres_freq <- get_thres_freq(cdf_data, threshold)
    thres_freq_lbl <- paste(round(thres_freq*100, 1), "%")
    thres_freq_lbl_inv <- paste(round(100 - thres_freq*100, 1), "%")


    max_y <- max(cdf_data$val)
    min_y <- min(cdf_data$val)
    print("thres freq:")
    print(thres_freq)

    #print(paste(expr_name, " median:", quantile(x)))
    #x <- x[which(x <= 5)]
    #hist(x, xlab=expr_name)
    #plot(ecdf(x), main=expr_name)
    #df <- cbind(x, y)
    fig <- ggplot(cdf_data, aes(x = val, y = freq))
    fig <- fig + geom_line()
    #fig <- fig + coord_cartesian(xlim = c(-0.5, quantile(x,c(0.90))[1]))
    fig <- fig + scale_x_log10()
    fig <- fig + ylab("Cumulative Frequency")
    fig <- fig + xlab("Expression Value")
    #fig <- fig + labs(title=expr_name)

    # TODO: change lty
    fig <- fig + annotate("segment", y = thres_freq, yend = thres_freq, x = min_y, xend = max_y, colour = "blue", lty="dotted")
    fig <- fig + annotate("segment", y = 0, yend = 1.0, x = threshold, xend = log(threshold,base=10), colour = "blue",lty="dotted")
    # text annotation
    fig <- fig + annotate("text", y=0, x=threshold, label=paste("t <",t_lbl), vjust=0, size=4)
    #fig <- fig + annotate("text", x=0, y=threshold, label=paste("t >=",threshold), hjust=0, vjust=-0.2, size=4)

    # text annotation to freq
    fig <- fig + annotate("text", y=thres_freq, x=min_y, label=thres_freq_lbl, vjust=1.1, hjust=0, size=4)
    # annotation lines


    return(fig)
}


# plots the percentage of expressed classifications linearly on the x-axis
# TODO: this should be called ` plot HPA expression diagram`
plot_expr_inv_cdf_nolog <- function(expr_name="hpa", threshold = 1, t_lbl=threshold)
{
    # get the data
    data <- get_expr_values(expr_name)
    
    #y <- unique(data$ExpressionValue)
    #
    #inv_cdf <- function(val)
    #{
    #    return (sum(data$ExpressionValue <= val) / length(data$ExpressionValue))
    #}

    # insert a (0,0) column at the beginning
    cdf_data <- inverse_cdf(data$ExpressionValue)
    #cdf_data <- rbind(cdf_data, c(0,0))

    # get the freq of the cdf at the threshold position
    thres_freq <- get_thres_freq(cdf_data, threshold)
    thres_freq_lbl <- paste(round(thres_freq*100, 1), "%")
    thres_freq_lbl_inv <- paste(round(100 - thres_freq*100, 1), "%")


    max_y <- max(cdf_data$val)
    min_y <- min(cdf_data$val)
    print("thres freq:")
    print(thres_freq)

    #print(paste(expr_name, " median:", quantile(x)))
    #x <- x[which(x <= 5)]
    #hist(x, xlab=expr_name)
    #plot(ecdf(x), main=expr_name)
    #df <- cbind(x, y)
    fig <- ggplot(cdf_data, aes(x = val, y = freq))
    fig <- fig + geom_line()
    fig <- fig + coord_cartesian(xlim = c(-0.5, max_y))
    fig <- fig + ylab("Cumulative Frequency")
    fig <- fig + xlab("Expression Value")
    #fig <- fig + labs(title=expr_name)

    # annotation lines
    # TODO: change lty
    fig <- fig + annotate("segment", y = thres_freq, yend = thres_freq, x = min_y, xend = max_y, colour = "blue", lty="dotted")
    fig <- fig + annotate("segment", y = 0, yend = 1.0, x = threshold, xend = threshold, colour = "blue",lty="dotted")
    # text annotation
    fig <- fig + annotate("text", y=0, x=threshold, label=paste("t <",t_lbl), vjust=0, size=4)
    #fig <- fig + annotate("text", x=0, y=threshold, label=paste("t >=",threshold), hjust=0, vjust=-0.2, size=4)

    # text annotation to freq
    fig <- fig + annotate("text", y=thres_freq, x=min_y, label=thres_freq_lbl, vjust=1.1, hjust=0, size=4)
    #fig <- fig + annotate("text", x=thres_freq, y=1, label=thres_freq_lbl_inv, vjust=0, hjust=0, size=4)


    return(fig)
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

save_plots_tissue_expr_distr <- function()
{
    pdf("../figs/tissue_expression_distr.pdf", width=6, height=4.6)
    p <- plot_all_expr(plot_tissue_expr_count_hist)
    print(p)
    dev.off()
}


save_plots_tshk_thresholds <- function()
{
    pdf("../figs/tshk_thresholds.pdf", width=6, height=3.5)
    p <- plot_all_expr(plot_tshk_thresholds)
    print(p)
    dev.off()
}

save_plots_tshk_example <- function()
{
    figs <- list()
    p1 <- plot_tshk_example("hpa", 0.15)
    figs <- c(figs, list(p1))

    p <- plot_tshk_example("hpa", 0.25)
    figs <- c(figs, list(p))



    par(mfrow=c(1,2))
    #all_figs <- do.call(grid.arrange, figs)
    #all_figs <- do.call(arrangeGrob, figs)

    # TODO: FIXME this still plots one plot rather than two
    pdf("../figs/tshk_class_example.pdf", width=6, height=2.6)
    all_figs <- do.call(grid.arrange, c(figs,list(nrow=1)))
    print(all_figs)
    dev.off()
    return (figs)
}

save_plots_tshk_all <- function()
{
    figs <- list()


    # Illumina body map:
    p <- plot_tshk_example("emtab", 0.15)
    p <- p + labs(title="Illumina Body Map 2.0")
    figs <- c(figs, list(p))

    # GeneAtlas
    p <- plot_tshk_example("gene_atlas", 0.15)
    p <- p + labs(title="Gene Atlas")
    figs <- c(figs, list(p))

    # RNAseq atlas
    p <- plot_tshk_example("rnaseq_atlas", 0.15)
    p <- p + labs(title="RNAseq Atlas")
    figs <- c(figs, list(p))

    # TODO: HPA
    p <- plot_tshk_example("hpa", 0.15)
    p <- p + labs(title="Human Protein Atlas")
    figs <- c(figs, list(p))


    par(mfrow=c(2,2))


    #all_figs <- do.call(grid.arrange, figs)
    all_figs <- do.call(arrangeGrob, figs)

    # actually plot:

    pdf("../figs/tshk_all.pdf", width=6, height=4.8)
    print(all_figs)
    dev.off()
}

save_plots_expr_value_hist <- function()
{
    figs <- list()

    # TODO: HPA
    p <- plot_expr_inv_cdf_nolog("hpa", 1.0)
    p <- p + labs(title="Human Protein Atlas")
    figs <- c(figs, list(p))

    # GeneAtlas
    p <- plot_expr_inv_cdf("gene_atlas", 100.0, "100")
    p <- p + labs(title="Gene Atlas")
    figs <- c(figs, list(p))

    # Illumina body map:
    p <- plot_expr_inv_cdf("emtab", 1.0, "1.0")
    p <- p + labs(title="Illumina Body Map 2.0")
    figs <- c(figs, list(p))

    # RNAseq atlas
    p <- plot_expr_inv_cdf("rnaseq_atlas", 1.0, "1.0")
    p <- p + labs(title="RNAseq Atlas")
    figs <- c(figs, list(p))


    par(mfrow=c(2,2))


    #all_figs <- do.call(grid.arrange, figs)
    all_figs <- do.call(arrangeGrob, figs)

    # actually plot:

    pdf("../figs/expression_thresholding.pdf", width=6, height=4.8)
    print(all_figs)
    dev.off()
}
#plot(all_figs)



