library(ggplot2)
library(gridExtra) # for `grid.arrange`

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
#print(all_figs)
print(all_figs)
dev.off()
#plot(all_figs)

