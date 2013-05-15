# TODO: Add comment
# 
# Author: flick
###############################################################################



scatterplot_by_degree <- function(full_ppi_graph, low_spec_properties, high_spec_properties, property_name)
{

    degrees <- degree(full_ppi_graph)

    # non_nan_knn <- graph_knn$knn[!is.na(graph_knn$knn)]
    degrees_high <- degrees[match(names(high_spec_properties), names(degrees))]
    degrees_low  <- degrees[match(names(low_spec_properties), names(degrees))]

    # plot as dotplot 
    plot(degrees_high, high_spec_properties, type="p", col="red", xlab="Degree", ylab=property_name) #, log="xy")
    lines(degrees_low, low_spec_properties, type="p", col="blue")

    # TODO stddev of ALL data, meaning for every vertex, the degrees of all of his neighboors into a table
    # TODO only then aggregate with mean and std-dev
    # TODO especially no multi-stage mean (messes up the total std-dev)

    # draw legend
    legend("topright",legend=c("Specific Proteins", "Promiscuous Proteins"), fill=c("red","blue"))

}


mean_log_lm_plot_by_degree <- function(full_ppi_graph, low_spec_properties, high_spec_properties, property_name)
{

    degrees <- degree(full_ppi_graph)

    # throw away nodes with degree 0 (log plot)
    degrees = degrees[which(degrees != 0)]


    # non_nan_knn <- graph_knn$knn[!is.na(graph_knn$knn)]
    degrees_high <- degrees[match(names(high_spec_properties), names(degrees))]
    degrees_low  <- degrees[match(names(low_spec_properties), names(degrees))]

    # aggregate data
    high_spec_df <- as.data.frame(cbind(degrees_high, high_spec_properties))
    high_spec_means <- aggregate(high_spec_df["high_spec_properties"], by = high_spec_df["degrees_high"], FUN=mean)
    low_spec_df <- as.data.frame(cbind(degrees_low, low_spec_properties))
    low_spec_means <- aggregate(low_spec_df["low_spec_properties"], by = low_spec_df["degrees_low"], FUN=mean)

    # throw out those properties that have value 0 (log plot)
    high_spec_means <- high_spec_means[which(high_spec_means$high_spec_properties != 0),]
    low_spec_means <- low_spec_means[which(low_spec_means$low_spec_properties != 0),]



    out1 <<- high_spec_means
    out2 <<- low_spec_means

    plot(high_spec_means$degrees_high, high_spec_means$high_spec_properties, type="p", log="xy", col="red", xlab= "Degree", ylab=property_name)
    abline(lm(log10(high_spec_means$high_spec_properties)~log10(high_spec_means$degrees_high)),col="red", lty=2)

    lines(low_spec_means$degrees_low, low_spec_means$low_spec_properties, type="p", col="blue")
    abline(lm(log10(low_spec_means$low_spec_properties)~log10(low_spec_means$degrees_low)),col="blue", lty=2)

}



boxplot_by_degree <- function(full_ppi_graph, promiscuous_properties, specific_properties, property_name)
{
    
    degrees <- degree(full_ppi_graph)

    # non_nan_knn <- graph_knn$knn[!is.na(graph_knn$knn)]
    degrees_high <- degrees[match(names(specific_properties), names(degrees))]
    degrees_low  <- degrees[match(names(promiscuous_properties), names(degrees))]

    xaxis_lims = range(degrees_low, degrees_high)
    xaxis_lims[2] = 20

    boxplot(promiscuous_properties~degrees_low, at = sort(unique(degrees_low))-0.15, boxwex=0.3, xlim = xaxis_lims ,border="blue", xaxt="n", xlab="Degree", ylab=property_name)
    boxplot(specific_properties~degrees_high, at = sort(unique(degrees_high))+0.15, boxwex=0.3, border="red", add=TRUE, xaxt = "n")
    axis(1)
    legend("topright",legend=c("Specific Proteins", "Promiscuous Proteins"), fill=c("red","blue"))
}

