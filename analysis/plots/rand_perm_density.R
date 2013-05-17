


source("../rand_perm_tests.R", chdir=TRUE)

# does random permutation testing on the combined c(class1, class2), plots the densities
# and the actual values for class1 and class2 based on the value returned by the property_function
plot_rand_perm_density_2class <- function(class1, class2, property_function, num_perm_tests, property_name)
{
    # get all items to randomly select from
    all_items <- c(class1, class2)

    # get the property value for both classes
    property_class1 <- property_function(class1)
    property_class2 <- property_function(class2)

    # get the threshold value
    threshold <- length(class1)/(length(all_items))

    # get the density for random permutation tests for the property
    density_class1 <- rand_perm_test_density(all_items, property_function, num_perm_tests, threshold)
    density_class2 <- rand_perm_test_density(all_items, property_function, num_perm_tests, 1-threshold)

    # calc limits for the plot
    xlims <- range(density_class1$x,density_class2$x, property_class1, property_class2)
    ylims <- range(density_class1$y,density_class2$y)


    # plot the density of the assortativity of the random sampling
    plot(NaN, xlim=xlims, ylim=ylims, xlab=property_name, ylab="density of random samples")
    lines(density_class1, col="red")
    lines(density_class2, col="blue")
    abline(v = property_class1, col="red")
    abline(v = property_class2, col="blue")

}

