# returns all expresison data sets aliases
get_exprs <- function()
{
    exprs <- c("emtab", "gene_atlas", "rnaseq_atlas", "hpa") #, "hpa_all")
    return(exprs)
}

# returns short name representations of the expression data sets
to_short_expr_name <- function(s)
{
    if (length(s) > 1)
    {
        # if a vector/lsit/matrix is passed in, apply
        # this function recursively on all the elements
        return(sapply(s, to_short_expr_name))
    }

    if (s == "emtab")
    {
        return ("Body Map")
    }
    else if (s == "gene_atlas")
    {
        return ("Gene Atlas")
    }
    else if (s == "rnaseq_atlas")
    {
        return ("RNAseq Atlas")
    }
    else if (s == "hpa")
    {
        return ("HPA")
    }
    else if (s == "hpa_all")
    {
        return ("HPA All")
    }
    else
    {
        return ("ERROR ERROR ERROR FIXME")
    }
}

to_short_expr_name_v <- function(v)
{

}


# returns name representations of the expression data sets
to_expr_name <- function(s)
{
    if (length(s) > 1)
    {
        # if a vector/lsit/matrix is passed in, apply
        # this function recursively on all the elements
        return(sapply(s, to_expr_name))
    }

    if (s == "emtab")
    {
        return ("Illumina Body Map 2.0")
    }
    else if (s == "gene_atlas")
    {
        return ("Gene Atlas")
    }
    else if (s == "rnaseq_atlas")
    {
        return ("RNAseq Atlas")
    }
    else if (s == "hpa")
    {
        return ("Human Protein Atlas")
    }
    else if (s == "hpa_all")
    {
        return ("Human Protein Atlas (all)")
    }
    else
    {
        return ("ERROR ERROR ERROR FIXME")
    }
}

