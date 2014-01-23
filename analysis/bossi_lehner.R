
# load required libs:
library(ggplot2)
library(plyr) # for `round_any`
library(gridExtra) # for `grid.arrange`
library(Hmisc) # for spearman correlation + p-value

# returns a list of all PPIs that can be analyzed
# TODO: put this as a unified function somewhere (using the SQL database)
# TODO: also so that the order of PPIs is always the same for all plots
get_ppis <- function()
{
    # FIXME: do this properly (not hardcoded)
    ppis <- c("bossi", "string", "psicquic_all", "havu", "ccsb")
    return(ppis)
}

get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa", "hpa_all")
    return(exprs)
}

###########################################
#  functions for calculating error bars:  #
###########################################


# standard error:
std <- function(x) sd(x)/sqrt(length(x))

# mean plus/minus one standard error
upper_se <- function(x) mean(x) + std(x)
lower_se <- function(x) mean(x) - std(x)



##########################
#  Bossi & Lehner Fig1B  #
##########################


get_expression_property_data <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="degree_test")
{

    # load sql config and get connection
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # create table names
    prop_table <- paste("prop", ppi_name, expr_name, prop_name, sep="_")
    expr_count_table <- paste(expr_name,"core_expr_counts", sep="_")

    # create SQL query to get ExpressionCount vs Degree
    query <- paste("SELECT a.Gene, a.Value, b.ExpressedCount, b.TotalCount
                    FROM ", prop_table, " AS a",
                   "INNER JOIN ", expr_count_table, " AS b",
                   "ON a.Gene = b.Gene")

    # load ppi network from db
    data <- dbGetQuery(con, query)
    return(data);
}


# reproducing the Bossi&Lehner plot (Fig. 1 B)
plot_degree_vs_expr <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="maxts_degree")
{

    data <- get_expression_property_data(ppi_name, expr_name, prop_name)

    # round off the exressedcount as simple "binning"
    # TODO: might be able to do this better :)
    max_expr_count <- max(data$ExpressedCount)
    round_to <- 10
    if (max_expr_count < 60)
    {
        round_to <- round(max_expr_count / 7)
    }
    #data$ExpressedCount <- round_any(data$ExpressedCount, round_to)

    error_width <- max_expr_count / 20

    print(paste(ppi_name, expr_name,":",length(data$ExpressedCount)))

    # plot using ggplot2
    fig <- ggplot(data, aes(x=ExpressedCount, y=Value)) +
          stat_summary(fun.y="mean", geom="line") +
          stat_summary(fun.y="mean", geom="point") +
          # draw error bars:
          stat_summary(fun.ymin="lower_se", fun.ymax="upper_se",
                       geom="errorbar", width=error_width) +
          labs(title=paste("PPI:", ppi_name, ", Expr:", expr_name),
               xlab="Number of tissues in which protein is expressed",
               ylab="Protein interaction degree")

    return(fig)
}


# checking for bossi results:
# 1.) tissue specific proteins make fewer protein interactions than widely
#     expressed proteins (use spearman)
# 2.) most recently evolved proteins have fewer interactions than ancient proteins TODO
test_bossi_1 <- function(ppi_name="bossi", expr_name="gene_atlas", prop_name="maxts_degree")
{
    data <- get_expression_property_data(ppi_name, expr_name, prop_name)

    # run a spearman test:

    for (test in c("spearman", "kendall", "pearson"))
    {
        tst = cor.test(data$ExpressedCount, data$Value, method=test)
        rho = tst$estimate
        pvalue = tst$p.value
        print(paste("runing", test, "test: ", ppi_name, expr_name, "rho=", rho, ", p-value=", pvalue))
    }
}

# 3.) to which extend to tissue-specific proteins interact with most widely expressed proteins?
test_bossi_3 <- function(ppi_name="bossi", expr_name="gene_atlas", prop_name="maxts_degree")
{
    # needs data from the actual graph -> TODO implement in NetworKit!?
}


create_test_degrees <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="degree_test")
{
    
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # create table names
    prop_table <- paste("prop", ppi_name, expr_name, prop_name, sep="_")
    
    dbSendQuery(con, paste("DROP TABLE IF EXISTS ", prop_table))
    query <- paste("CREATE TABLE IF NOT EXISTS ", prop_table, " AS ",
                   "SELECT Gene, COUNT(*) AS Value FROM ",
                   "(",
                   "SELECT Gene1 AS Gene FROM ", ppi_name,
                   # requires UNION all, otherwise all degrees will be 1
                   " UNION ALL ",
                   "SELECT Gene2 As Gene FROM ", ppi_name,
                   ")",
                   "GROUP BY Gene")
    dbSendQuery(con, query);
}


compare_test_degrees <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="degree_test")
{
    
    source("sql_config.R")
    con <- get_sql_conn('/home/patrick/dev/bio/data/test_matching.sqlite')

    # create table names
    prop_table_1 <- paste("prop", ppi_name, expr_name, "_degree", sep="_")
    prop_table_2 <- paste("prop", ppi_name, expr_name, "degree_test", sep="_")

    query <- paste("SELECT a.Gene, a.Value AS Value1, b.Value AS Value2 ",
                   "FROM ", prop_table_1, " AS a INNER JOIN ",
                   prop_table_2, " AS b ON a.Gene = b.Gene")
    data <- dbGetQuery(con, query)

    #print(which(data$Value1 != data$Value2))
}

test_for_all <- function(func)
{
    ppis <- get_ppis()
    exprs <- get_exprs()
    for (p in ppis)
    {
        for (e in exprs)
        {
            func(p,e)
        }
    }
}


plot_all <- function()
{
    ppis <- get_ppis()
    exprs <- get_exprs()

    par(mfrow=c(length(ppis),length(exprs)))

    figs <- list()
    for (p in ppis)
    {
        for (e in exprs)
        {
            # call the bossi plotting function
            #fig <- plot_degree_vs_expr(p, e, "coexpr_degree")
            create_test_degrees(p,e);
            fig <- plot_degree_vs_expr(p, e, "degree_test")
            figs <- c(figs, list(fig))
        }
    }

    all_figs <- do.call(grid.arrange, figs)

    # in order to pass named parameters to grid.arrange:
    #do.call(grid.arrange, c(figs,list(nrow=5,left="hello there"))
}




