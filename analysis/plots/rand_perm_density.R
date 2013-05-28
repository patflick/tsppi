


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
    rand_class1 <- rand_perm_test(all_items, property_function, num_perm_tests, threshold)
    rand_class2 <- rand_perm_test(all_items, property_function, num_perm_tests, 1-threshold)

    density_class1 <- density(rand_class1)
    density_class2 <- density(rand_class2)

    # calc limits for the plot
    xlims <- range(density_class1$x,density_class2$x, property_class1, property_class2)
    ylims <- range(density_class1$y,density_class2$y)


    # plot the density of the assortativity of the random sampling
    plot(NaN, xlim=xlims, ylim=ylims, xlab=property_name, ylab="density of random samples")
    lines(density_class1, col="red")
    lines(density_class2, col="blue")
    abline(v = property_class1, col="red")
    abline(v = property_class2, col="blue")


    # calculate p values
    sd_class1 = sd(rand_class1)
    mean_class1 = mean(rand_class1)
    #pvalue_class1 = 1-pnorm(property_class1, mean=mean_class1, sd=sd_class1)
    pvalue_class1 = round(1-pnorm(abs(property_class1-mean_class1)/sd_class1), digits=4)
    pvalue_class1_text = paste(c("p ~= ", pvalue_class1), sep="", collapse="")
    # print(pvalue_class1_text)
    text(property_class1, 0, pvalue_class1_text, pos=4, col="red")

    # calculate p values
    sd_class2 = sd(rand_class2)
    mean_class2 = mean(rand_class2)
    #pvalue_class2 = 2-pnorm(property_class2, mean=mean_class2, sd=sd_class2)
    pvalue_class2 = round(1-pnorm(abs(property_class2-mean_class2)/sd_class2), digits=4)
    pvalue_class2_text = paste(c("p ~= ", pvalue_class2), sep="", collapse="")
    # print(pvalue_class2_text)
    text(property_class2, 0, pvalue_class2_text, pos=2, col="blue")
}

