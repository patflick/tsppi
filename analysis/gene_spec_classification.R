# loads the genes from the database in a sorted order and splits them
# according to a given threshold into high and low specific genes
# 
# Author: Patrick Flick
###############################################################################


# load utils
source("utils.R", chdir=TRUE)

# holds the HPA expression data as soon as it is loaded
# making sure the data is only loaded once
EXPRESSED_GENES_DATA <- 0

# loads the expressed genes data from the database in
# sorted order. then returns the data. if the data
# has already been loaded, it is not loaded from the
# database again. the data in memory is returned.
load_expression_ordered_genes <- function()
{
    if (EXPRESSED_GENES_DATA == 0)
    {
        # load sql config and get connection
        source("sql_config.R", chdir=TRUE)
        con <- get_sql_conn()
    
        EXPRESSED_GENES_DATA <<- dbGetQuery(con, "
        SELECT
            Gene,
            CountHigh,
            CountMedium,
            CountLow,
            CountExpressed,
            CountTotal
        FROM hpa_gene_levels
        ORDER BY CountHigh+CountMedium+CountLow ASC
        ")
    
        dbDisconnect(con)
    }

    return(EXPRESSED_GENES_DATA)
}


# load and returns the genes of the HPA in two classes as list[promiscuous_genes, specific_genes]
get_genes_in_specificity_classes <- function(threshold=0.5)
{
    # get the gene data
    data <- load_expression_ordered_genes()

    # get number of genes
    n <- nrow(data)

    # get the high and low specificity genes
    specific_genes <- data$Gene[1:floor(threshold*n)]
    promiscuous_genes <- data$Gene[(floor(threshold*n)+1):n]

    # return a list with the two gene vectors (promiscuous and specific)
    return(list(promiscuous_genes=promiscuous_genes, specific_genes=specific_genes))
}


# returns the given vertex properties of the network, classified into two
# gene specificity classes as list[promiscuous_properties, specific_properties]
# according to the classification threshold given.
get_vertex_properties_in_spec_classes <- function(vertex_property, threshold = 0.5)
{
    # get the genes
    list[promiscuous_genes, specific_genes] <- get_genes_in_specificity_classes(threshold)

    # get the low and high spec properties from the selected genes
    promiscuous_properties = vertex_property[which(names(vertex_property) %in% promiscuous_genes)]
    specific_properties = vertex_property[which(names(vertex_property) %in% specific_genes)]
    
    # return a list with the two gene vectors (promiscuous and specific)
    return(list(promiscuous_properties=promiscuous_properties, specific_properties=specific_properties))
}



