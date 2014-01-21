
# load required libs:
library(ggplot2)
library(plyr) # for `round_any`
library(gridExtra) # for `grid.arrange`

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


# reproducing the Bossi&Lehner plot (Fig. 1 B)
plot_degree_vs_expr <- function(ppi_name="bossi",expr_name="gene_atlas",prop_name="maxts_degree")
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

    # round off the exressedcount as simple "binning"
    # TODO: might be able to do this better :)
    max_expr_count <- max(data$ExpressedCount)
    round_to <- 10
    if (max_expr_count < 60)
    {
        round_to <- round(max_expr_count / 7)
    }
    data$ExpressedCount <- round_any(data$ExpressedCount, round_to)

    error_width <- max_expr_count / 20

    # standard error:
    std <- function(x) sd(x)/sqrt(length(x))

    # mean plus/minus one standard error
    upper_se <- function(x) mean(x) + std(x)
    lower_se <- function(x) mean(x) - std(x)

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


ppis <- get_ppis()
exprs <- get_exprs()

par(mfrow=c(length(ppis),length(exprs)))

figs <- list()
for (p in ppis)
{
    for (e in exprs)
    {
        # call the bossi plotting function
        fig <- plot_degree_vs_expr(p, e, "coexpr_degree")
        figs <- c(figs, list(fig))
    }
}

do.call(grid.arrange, figs)
#grid.arrange(figs, ncol=length(exprs))


