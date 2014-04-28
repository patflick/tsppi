

library(reshape) # for `melt`A
library(ggplot2)
library(xtable) # for latex tables

source("./ppi_utils.R")
source("./expr_utils.R")

get_cc_data <- function () {
    data <- read.csv("./data/benchmark_cc.csv", header=TRUE, sep=";")

    # sort by naive time
    data <- data[with(data, order(-naive)),]

    # map to short name
    data$ppi <- to_short_ppi_name(data$ppi)
    data$expr <- to_short_expr_name(data$expr)

    data$pe <- paste(data$ppi, "-", data$expr)

    return(data)
}

acc_cc_data <- function()
{
    data <- get_cc_data();

    # create empty data frame
    res <- data.frame(Method=c("NetworKit", "Neighbor combinations", "Tissue expr. vectors"), Full.Runtime=rep(0.0, 3))

    res$Full.Runtime[1] = round(sum(data$naive),1)
    res$Full.Runtime[2] = round(sum(data$neighcomb),1)
    res$Full.Runtime[3] = round(sum(data$tsvector),1)

    return(res)
}


plot_cc_benchmark <- function()
{
    data <- get_cc_data()

    plot_data <- data[1:5,c("pe", "naive", "neighcomb", "tsvector")]
    colnames(plot_data) <- c("PPI.Expr", "NetworKit", "Neighbor combinations", "Tissue expr. vectors")

    # for plotting
    plot_data <- melt(plot_data, id=c("PPI.Expr"))
    colnames(plot_data) <- c("PPI.Expr", "Method", "Runtime")

    fig <- ggplot(plot_data, aes(PPI.Expr, Runtime, fill=Method)) +
            geom_bar(stat="identity", position="dodge") +
            geom_text(aes(PPI.Expr, pmin(20,Runtime),  label=paste(round(Runtime,2),"s"), hjust=0), size=3.5, position = position_dodge(width=1)) +
            coord_flip(ylim=c(0,25)) +
            xlab("") +
            ylab("Run time [s]") +
            labs(title="Run time of methods for computation of clustering coefficients")

    return (fig)
}

get_bw_data <- function(benchmark="omp_dynamic")
{
    # get serial benchmark
    data <- read.csv(paste("./data/benchmark_bw_",benchmark ,".csv", sep=""), header=FALSE, sep=";")
    colnames(data) <- c("ppi", "expr", "outerloop", "innerloop")

    # sort by naive time
    data <- data[with(data, order(-outerloop)),]

    # map to short name
    data$ppi <- to_short_ppi_name(data$ppi)
    data$expr <- to_short_expr_name(data$expr)

    data$pe <- paste(data$ppi, "-", data$expr)

    return (data)
}

acc_bw_data <- function()
{
    par_data <- get_bw_data("omp_dynamic")
    seq_data <- get_bw_data("serial")

    res <- data.frame(Method=c("Create Subgraphs", "Use Tissue Vectors"),
                      Sequential=rep(0.0,2), Parallel=rep(0.0,2),
                      Speedup=rep(0.0,2))

    res$Sequential[1] = sum(seq_data$outerloop)
    res$Sequential[2] = sum(seq_data$innerloop)

    res$Parallel[1] = sum(par_data$outerloop)
    res$Parallel[2] = sum(par_data$innerloop)

    res$Speedup = round(res$Sequential / res$Parallel, 2)

    res$Sequential = round(res$Sequential,1)
    res$Parallel = round(res$Parallel,1)

    return(res)
}

plot_bw_benchmark <- function()
{
    data <- get_bw_data()

    plot_data <- data[1:8,c("pe", "outerloop", "innerloop")]
    colnames(plot_data) <- c("PPI.Expr", "Create Subgraphs", "Use Tissue Vectors")

    # for plotting
    plot_data <- melt(plot_data, id=c("PPI.Expr"))
    colnames(plot_data) <- c("PPI.Expr", "Method", "Runtime")

    fig <- ggplot(plot_data, aes(PPI.Expr, Runtime, fill=Method)) +
            geom_bar(stat="identity", position="dodge") +
            geom_text(aes(PPI.Expr, pmin(max(Runtime)*.60,Runtime),  label=paste(round(Runtime,1),"s"), hjust=0), size=3.5, position = position_dodge(width=1)) +
            coord_flip() +
            xlab("") +
            ylab("Run time [s]") +
            labs(title="Run time of methods for computation of clustering coefficients")
    return (fig)
}

save_plots <- function()
{
    # clustering coeff. benchmark
    fig <- plot_cc_benchmark()
    pdf("../figs/benchmark_cc.pdf", width=8, height=3.2)
    print(fig)
    dev.off()

    # betweenness benchmarks
    fig <- plot_bw_benchmark()
    pdf("../figs/benchmark_bw.pdf", width=8, height=3.5)
    print(fig)
    dev.off()
}
