
rand_vertex_subset_property <- function(all_vertices, subset_property_function, subset_size=0.5)
{
    n <- length(all_vertices)
    # creates a random permutation of the genes
    vertex_subset <- sample(all_vertices, floor(subset_size*n))


    # calculate property of random subset
    subset_property <- subset_property_function(vertex_subset)
    return(subset_property)
}

rand_perm_test_density <- function(all_vertices, subset_property_function, repeats=10000, subset_size=0.5)
{
    # do the random permutation test, and get the list of values
    values <- rand_perm_test(all_vertices, subset_property_function, repeats, subset_size)

    # return the density of the values
    return(density(na.omit(values)))
}


rand_perm_test <- function(all_vertices, subset_property_function, repeats=10000, subset_size=0.5)
{
    values <- c()
    for (i in 1:repeats)
    {
        # add new value to vector
        values <- c(values, rand_vertex_subset_property(all_vertices, subset_property_function, subset_size))

        # if part of 20% rise, update console
        if (i %% floor(repeats/5) == 0)
        {
            # TODO print the percentage
            print(i)
        }
    }
    # return the values (results) of the random tests
    return(values)
}


