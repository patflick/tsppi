

library(reshape) # for `melt`A
library(ggplot2)

source("./ppi_utils.R")
source("./expr_utils.R")


data <- read.csv("./data/benchmark_cc.csv", header=TRUE, sep=";")

# sort by naive time
data <- data[with(data, order(-naive)),]

# map to short name
data$ppi <- to_short_ppi_name(data$ppi)
data$expr <- to_short_expr_name(data$expr)

data$pe <- paste(data$ppi, "-", data$expr)

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
        #coord_cartesian(ylim = c(0, 25))

pdf("../figs/benchmark_cc.pdf", width=8, height=3.2)
print(fig)
dev.off()
