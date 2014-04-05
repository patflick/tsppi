
library(ggplot2)
library(xtable) # for converting a data frame into a latex table
library(igraph) # for power.law.fit

get_ppis <- function()
{
    # FIXME: do this properly (not hardcoded)
    ppis <- c("bossi", "string", "psicquic_all", "havu", "ccsb")
    return(ppis)
}

to_short_expr_name <- function(s)
{
    if (s == "bossi")
    {
        return ("Bossi")
    }
    else if (s == "string")
    {
        return ("STRING")
    }
    else if (s == "psicquic_all")
    {
        return ("IMEx")
    }
    else if (s == "havu")
    {
        return ("Havugimana")
    }
    else if (s == "ccsb")
    {
        return ("HI-2012")
    }
    else
    {
        return ("ERROR ERROR ERROR FIXME")
    }
}

to_expr_name <- function(s)
{
    if (s == "bossi")
    {
        return ("Bossi & Lehner")
    }
    else if (s == "string")
    {
        return ("STRING DB")
    }
    else if (s == "psicquic_all")
    {
        return ("IMEx PSICQUIC")
    }
    else if (s == "havu")
    {
        return ("Havugimana et al.")
    }
    else if (s == "ccsb")
    {
        return ("CCSB HI-2012")
    }
    else
    {
        return ("ERROR ERROR ERROR FIXME")
    }
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
    degrees <- node_props$degree
    props[["avg_degree"]] <- mean(degrees)
    props[["max_degree"]] <- max(degrees)

    local_cc <- node_props$ClusteringCoeff
    props[["avg_local_cc"]] <- mean(local_cc)


    # return result
    return(props)
}


get_ppi_global_stats_table <- function()
{
    stats_tbl <- data.frame()
    for (p in get_ppis())
    {
        stats <- get_ppi_graph_stats(p)
        stats$ppi <- p
        stats_tbl <- rbind(stats_tbl, t(stats))
    }
    return(stats_tbl)
}

get_latex_summary_stats <- function()
{
    # get data and map alias to actual PPI name
    data <- get_ppi_global_stats_table()
    data$ppi_name <- sapply(data$ppi, to_short_expr_name)

    # select columns to output
    out_data <- data[,c("ppi_name", "n", "m", "avg_deg", "max_deg", "conn_comp")]

    # actually output the latex table to be copy-pasted into the report
    xtable(out_data)
}

# TODO: degree distribution fit
#   - print fit of degree distr to binomial, poisson and scale-free
#   - and according distribution parameters
#   - maybe plot showing the different degree distributions and the real one
#     (to show the actual goodness of fit)

fit_binom_degree_distr <- function(degrees)
{
    n <- length(degrees)
    m <- sum(degrees) / 2
    k <- max(degrees)
    p <- 2 * n / (n*(n-1))
    max_degr = k

    # Random graph model with binomial distribution
    x <- 1:max_degr
    freq <- dbinom(x, n-1, p)
    binom_model <- data.frame(degree=x, Freq=freq)

    return(binom_model)
}

fit_poisson_degree_distr <- function(degrees)
{
    avg_deg <- mean(degrees)
    # the lambda of the poisson distribution is the average degree
    lambda <- avg_deg

    # Random graph model with binomial distribution
    x <- 1:max_degr
    freq <- dpois(x, lambda)
    pois_model <- data.frame(degree=x, Freq=freq)
    return(pois_model)
}

fit_powerlaw_degree_distr <- function(degrees)
{
    # use igraph power.law.fit
    fit <- power.law.fit(degrees)

    return(fit)
}

# plot degree distribution of the networks
plot_degree_distr <- function(ppi_name="string")
{
    # get the data
    data <- get_node_properties(ppi_name)

    # get the degrees
    degree <- as.integer(data$degree)

    degree_distr <- as.data.frame(prop.table(table(degree)))
    degree_distr$degree <- as.integer(degree_distr$degree)


    # random graph (Erdos Renyi, or Gilbert)

    # get parameters for Gilbert random graph model G(N,p)
    # get binomial random model
    binom_model <- fit_binom_degree_distr(degree)
    # cut of zero freq
    binom_model <- binom_model[which(binom_model$Freq > 1e-20),]

    # poisson model
    pois_model <- fit_poisson_degree_distr(degree)

    # cut of zero freq
    pois_model <- pois_model[which(pois_model$Freq > 1e-20),]


    fig <- ggplot(degree_distr, aes(x=degree, y=Freq)) + geom_point(colour="grey60")
    fig <- fig + scale_x_log10() + scale_y_log10()
    fig <- fig + geom_line(data=binom_model, aes(x=degree, y=Freq), lty="dashed")
    fig <- fig + geom_line(data=pois_model, aes(x=degree, y=Freq), lty="dotted")
    fig <- fig + coord_cartesian(ylim =c(min(degree_distr$Freq)/5, max(degree_distr$Freq) * 5))


    return(fig)
}




# TODO: plot clustering coeff against p (N-1) [ see barabasi book ]
