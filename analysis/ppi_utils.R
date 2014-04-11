
get_ppis <- function()
{
    # FIXME: do this properly (not hardcoded)
    ppis <- c("bossi", "string", "psicquic_all", "havu", "ccsb")
    return(ppis)
}

to_short_ppi_name <- function(s)
{
    if (length(s) > 1)
    {
        # if a vector/lsit/matrix is passed in, apply
        # this function recursively on all the elements
        return(sapply(s, to_short_ppi_name))
    }

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

to_ppi_name <- function(s)
{
    if (length(s) > 1)
    {
        # if a vector/lsit/matrix is passed in, apply
        # this function recursively on all the elements
        return(sapply(s, to_ppi_name))
    }
    if (s == "bossi")
    {
        return ("Bossi & Lehner")
    }
    else if (s == "string")
    {
        return ("STRING")
    }
    else if (s == "psicquic_all")
    {
        return ("IMEx (PSICQUIC)")
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


plot_all_ppis <- function(plot_func, addTitle=FALSE, ...)
{
    figs <- list()
    for (p in get_ppis())
    {
        fig <- plot_func(p, ...)
        if (addTitle)
        {
            fig <- fig + labs(title=to_short_ppi_name(p))
        }
        figs <- c(figs, list(fig))
    }

    figs <- c(figs, list(degr_distr_legend()))

    all_figs <- do.call(grid.arrange, figs)
    return(all_figs)
}
