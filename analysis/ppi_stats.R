
library(ggplot2)

get_ppis <- function()
{
    # FIXME: do this properly (not hardcoded)
    ppis <- c("bossi", "string", "psicquic_all", "havu", "ccsb")
    return(ppis)
}


get_node_properties <- function(ppi_name = "string")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    ppi_node_prop_tbl <- paste(ppi_name, "node_properties", sep="_")
    query <- paste("SELECT * FROM ", ppi_node_prop_tbl)
    node_props <- dbGetQuery(con, query)

    return(node_props)
}

get_graph_properties <- function(ppi_name="string")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    ppi_prop_tbl <- paste(ppi_name, "properties", sep="_")
    query <- paste("SELECT * FROM ", ppi_prop_tbl)
    props <- dbGetQuery(con, query)

    return (props)
}


# TODO: print out graph information:
#   - n
#   - m
#   - connected components
#   - maximum degree
#   - average degree
#   - global clustering coeff
#   - avg local clustering coeff
#   - (?) degree assotativity 
get_ppi_graph_stats <- function(ppi_name="string")
{
    # will get and/or calculate all graph properties and save them into
    # named list elements
    props = list()

    # get the graph information from the database
    graph_props <- get_graph_properties(ppi_name)
    node_props <- get_node_properties(ppi_name)

    # extract easy info
    for (i in 1:dim(graph_props)[1])
    {
        props[[graph_props[i,1]]] = graph_props[i,2]
    }

    # TODO: summarize node wise data


    # return result
    return(props)
}


ppi_stats_table <- function()
{
    for (p in get_ppis())
    {
        stats <- get_ppi_graph_stats(p)
        # todo output as table! (i.e. somehow reformat as data frame and then
        # output as LaTeX table with the columns that i want!
        print(stats)
    }
}


# TODO: degree distribution fit
#   - print fit of degree distr to binomial, poisson and scale-free
#   - and according distribution parameters
#   - maybe plot showing the different degree distributions and the real one
#     (to show the actual goodness of fit)


# TODO: plot clustering coeff against p (N-1) [ see barabasi book ]
