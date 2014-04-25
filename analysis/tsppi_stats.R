
source("ppi_utils.R")
source("expr_utils.R")

get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all")
    return(exprs)
}




get_graph_properties <- function(ppi_name="string", expr_name="gene_atlas")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    ppi_prop_tbl <- paste(ppi_name, expr_name, "properties", sep="_")
    query <- paste("SELECT * FROM ", ppi_prop_tbl)
    props <- dbGetQuery(con, query)

    return (props)
}

get_ts_graph_properties <- function()
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    ppi_prop_tbl <- "ts_graph_properties"
    query <- paste("SELECT ppi, expr, Property, AVG(Value) as avg_value FROM ",
                   ppi_prop_tbl, " GROUP BY ppi, expr, Property")
    props <- dbGetQuery(con, query)
    return (props)
}

get_all_ts_properties <- function()
{
    data <- get_ts_graph_properties()

    properties <- unique(data$Property)

    # get first columns
    df <- data[which(data$Property == properties[1]),c("ppi","expr", "avg_value")]
    colnames(df) <- c("ppi", "expr", paste("ts_", properties[1], sep=""))

    for (i in 2:length(properties))
    {
        # add all the other columns
        new_col <- data[which(data$Property == properties[i]),c("avg_value")]

        # old names
        names <- colnames(df)
        df$new_col <- new_col
        colnames(df) <- c(names, paste("ts_", properties[i], sep=""))
    }

    return(df)
}

get_all_properties <- function()
{
    df <- data.frame()
    for (p in get_ppis())
    {
        for (e in get_exprs()) {
            x <- get_graph_properties(p, e)
            props <- t(x$Value)
            colnames(props) <- t(x$Property)
            props <- as.data.frame(props)
            props$ppi <- p
            props$expr <- e
            df <- rbind(df, props)
        }
    }
    return (df)
}


get_merged_properties <- function()
{
    x <- get_all_properties()
    y <- get_all_ts_properties()
    z <- merge(x,y)
    return(z)
}

# TODO: plot and percentages into grid plot
