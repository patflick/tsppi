source("./expr_utils.R")

get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all")
    return(exprs)
}




is_core <- function(expr_name="hpa")
{
    # load the ts/hk summary data from the database
    source("sql_config.R")
    con <- get_sql_conn()

    query <- paste("SELECT COUNT(DISTINCT Gene) FROM ", expr_name)
    gene_count <- dbGetQuery(con, query)
    query <- paste("SELECT COUNT(DISTINCT Type) FROM ", expr_name)
    tissue_count <- dbGetQuery(con, query)

    query <- paste("SELECT Gene, COUNT(DISTINCT Type) as cnt FROM ",
                   expr_name, " GROUP BY Gene HAVING cnt != ",
                   tissue_count)
    incompl_genes <- dbGetQuery(con, query)

    query <- paste("SELECT Type, COUNT(DISTINCT Gene) as cnt FROM ",
                   expr_name, " GROUP BY Type HAVING cnt != ",
                   gene_count)
    incompl_tissues <- dbGetQuery(con, query)


    return(list(genes=incompl_genes, tissues=incompl_tissues))
}


all_is_core <- function()
{
    for (e in get_exprs())
    {
        incompl <- is_core(e)
        ngenes <- dim(incompl$genes)[1]
        ntissues <- dim(incompl$tissues)[1]

        if (ngenes == 0 && ntissues == 0)
        {
            print(paste(e, " is complete"))
        } else {
            print(paste(e, " is incomplete with genes=", ngenes,
                        ", tissues=", ntissues))
        }
    }
}
