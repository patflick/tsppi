library(ggplot2)

source("ppi_utils.R")
source("expr_utils.R")

# plots tiles with text inside and scale based fill colouring
# for all PPI and EXPR
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

    # translate the db aliases into the real names of PPIs and Expression sets
    plot_data$expr <- to_short_expr_name(plot_data$expr)
    plot_data$ppi <- to_short_ppi_name(plot_data$ppi)

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

