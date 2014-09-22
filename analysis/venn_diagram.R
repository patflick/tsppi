######################################################################
# Creates venn diagrams from data stored in the database
#
# Needs the package 'Vennerable', install with:
#   source("http://bioconductor.org/biocLite.R")
#   biocLite(c("graph", "RBGL", "gtools", "xtable"))
#   install.packages("Vennerable", repos="http://R-Forge.R-project.org")
######################################################################


plot_venn_diagram <- function(agg_table)
{

    # load sql config and get connection
    source("sql_config.R")
    con <- get_sql_conn()


    # load ppi network from db
    venn_data <- dbGetQuery(con, paste(c("
    SELECT *
    FROM ", agg_table), sep="", collapse=""))

    # disconnect from db
    dbDisconnect(con)

    ## get all columns to use in the venn diagram
    cols <- colnames(venn_data)
    cols <- cols[which(cols != "count")]

    #venn_combs <- 1:nrow(venn_data)
    #venn_counts <- 1:nrow(venn_data)

    #for (i in 1:nrow(venn_data))
    #{
    #    bin_row <- venn_data[i,cols] == 1
    #    venn_combs[i] <- paste(cols[which(bin_row)], sep="", collapse="&")
    #    venn_counts[i] <- venn_data[i,"count"]
    #}

    #print(venn_combs)
    #print(venn_counts)

    # needs the package venneuler
    #library(venneuler)

    #vd <- venneuler(venn_combs, venn_counts)
    #plot(vd)

    library(Vennerable)

    venn_weights <- c()
    for (i in 1:nrow(venn_data))
    {
        bin_index <- paste(venn_data[i, cols], sep="", collapse="")
        venn_weights[bin_index] <- venn_data[i, "count"]
    }

    venn_obj <- Venn(SetNames=cols, Weight=venn_weights)

    if (length(cols) > 3)
    {
        plot(venn_obj, doWeights=FALSE, type="circles")
    }
    else
    {
        plot(venn_obj)
    }

    return(venn_obj)
}

