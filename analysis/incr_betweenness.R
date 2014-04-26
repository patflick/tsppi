
source("ppi_utils.R")
source("expr_utils.R")


SQL_INCR_TABLE_NAME <- "incr_betweenness"


get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all")
    return(exprs)
}

# analysis that determines proteins that become more important (according to
# betweenness centrality) in tissue specific networks

get_incr_query <- function(ppi_name="ccsb", expr_name="hpa")
{
    # get all table names
    ts_node_prop_tbl <- paste(ppi_name, expr_name, "ts_node_properties", sep="_")
    node_prop_tbl <- paste(ppi_name, expr_name, "node_properties", sep="_")
    expr_tbl <- paste(expr_name, "core", sep="_")
    expr_counts_tbl <- paste(expr_name, "core_expr_counts", sep="_")

    # concatenate the query
    query <- paste("SELECT '", ppi_name , "' AS ppi, '", expr_name, "' AS expr,
                   a.Gene, a.Tissue, a.Betweenness AS ts_betweenness,
                   c.Betweenness as betweenness, a.degree AS ts_degree,
                   c.degree AS degree, d.ExpressedCount, d.TotalCount
            FROM ", ts_node_prop_tbl, " AS a
            INNER JOIN ", expr_tbl, " as B ON a.Gene =  b.Gene AND a.Tissue = b.Type
            INNER JOIN ", node_prop_tbl, " AS c ON a.Gene = c.Gene
            INNER JOIN ", expr_counts_tbl, " AS d ON a.Gene = d.Gene
        WHERE b.Expressed = 1 AND a.Betweenness > c.Betweenness", sep="")

    return (query)
}

create_union_incr_table <- function()
{
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    queries <- c()
    for (p in get_ppis())
    {
        for (e in get_exprs())
        {
            # TODO
            q <- get_incr_query(p, e)
            queries <- c(queries, q)
        }
    }

    # concatenate all by union
    query <- paste(queries, collapse="  UNION  ")

    # delete old table
    dbSendQuery(con, paste("DROP TABLE IF EXISTS ", SQL_INCR_TABLE_NAME))
    dbSendQuery(con, paste("CREATE TABLE ", SQL_INCR_TABLE_NAME, " AS ", query))
}

get_incr_betweenness_nodes <- function()
{
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    query <- paste("SELECT * FROM ", SQL_INCR_TABLE_NAME)
    data <- dbGetQuery(con, query)

    return(data)
}

plot_betw_incr_distr <- function()
{
    data <- get_incr_betweenness_nodes()
}

# TODO:
# - further aggregation of data for analysis
# - find genes that "profit" in tissue specific networks
# - highest scoring: -> get some biological background on these (maybe they are good genes)
# - obviously this method has its limitations, since the betweenness centrality is also close
#   to being similar to a "bottleneck" in the network
