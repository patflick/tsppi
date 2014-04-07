
library(ggplot2)
library(xtable) # for converting a data frame into a latex table
library(gridExtra) # for `grid.arrange`
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
            fig <- fig + labs(title=to_short_expr_name(p))
        }
        figs <- c(figs, list(fig))
    }

    figs <- c(figs, list(degr_distr_legend()))

    all_figs <- do.call(grid.arrange, figs)
    return(all_figs)
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

dpowerlaw <- function(x, alpha, xmin=1)
{
    pdf <- ((alpha - 1)/xmin) * (x / xmin)**(-alpha)
    return (pdf)
}

ppowerlaw <- function(x, alpha, xmin=1)
{
    # integral over pdf starting from xmin
    # http://www.wolframalpha.com/input/?i=integral+from+b+to+t+of+%28%28a-1%29%2Fb%29%28x%2Fb%29%5E%28-a%29
    cdf <- 1 - (x / xmin)**(-alpha+1)
    return(cdf)
}


#######################################################################
#           Fit degree distributions (parameter estimation)           #
#######################################################################

fit_binom <- function(degrees)
{
    n <- length(degrees)
    m <- sum(degrees) / 2
    p <- 2 * m / (n*(n-1))

    params <- list(size=n-1, prob=p)
    return(params)
}

fit_poisson <- function(degrees)
{
    avg_deg <- mean(degrees)
    # the lambda of the poisson distribution is the average degree
    lambda <- avg_deg
    params <- list(lambda=lambda)
    return(params)
}

fit_exp <- function(degrees)
{
    rate <- 1 / mean(degrees)
    params <- list(rate=rate)
    return(params)
}

fit_powerlaw <- function(degrees)
{
    if (TRUE)
    {
        fit <- power.law.fit(degrees, implementation="R.mle")
        alpha <- coef(fit)
        xmin <- min(degrees)
    }
    else
    {
        fit <- power.law.fit(degrees, implementation="plfit")
        alpha <- fit$alpha
        xmin <- fit$xmin
    }
    params <- list(alpha=alpha, xmin=xmin)
    return(params)
}


#######################################################################
#              Fit degree distribution to 1:max(degrees)              #
#######################################################################

# fits the given degree distribution with the given fitter
# and applies the pdf with the fitted parameters to the range 1:max(degrees)
# and returns these values
fit_degree_distr <- function(degrees, fitter, pdf)
{
    # get pdf parameters
    params <- fitter(degrees)
    # get maximum degree and with it the input range
    max_degr <- max(degrees)
    x <- 1:max_degr

    # get the pdf at the given points
    freq <- do.call(pdf, c(list(x), params))
    model <- data.frame(degree=x, Freq=freq)

    return(model)
}


#######################################################################
#                  wrapper functions for convenience                  #
#######################################################################

fit_binom_degree_distr <- function(degrees)
{
    return(fit_degree_distr(degrees, fit_binom, dbinom))
}


fit_poisson_degree_distr <- function(degrees)
{
    return(fit_degree_distr(degrees, fit_poisson, dpois))
}

fit_powerlaw_degree_distr <- function(degrees)
{
    return(fit_degree_distr(degrees, fit_powerlaw, dpowerlaw))
}

fit_exp_degree_distr <- function(degrees)
{
    return(fit_degree_distr(degrees, fit_exp, dexp))
}


#######################################################################
#                      Plotting: generate legend                      #
#######################################################################

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

degr_distr_legend <- function()
{
    # TODO: create some plot using all the previous data
    x <- 1:10
    df <- data.frame(x=x, y=x, z=as.character(x), y1=x, y2=x, y3=x)

    dists <- c("Actual Distribution", "Binomial Distribution",
               "Exponential Distribution", "Power-law Distribution")
    x11()
    fig <- ggplot(data=df, aes(x=x, y=y, lty="a", color="a", shape="a")) + geom_point()+
            geom_line(data=df, aes(x=x,y=-y1/2,lty="b", color="b", shape="b")) +
            geom_line(data=df, aes(x=x,y=-y2,lty="c", color="c", shape="c")) +
            geom_line(data=df, aes(x=x,y=y3,lty="d", color="d", shape="d")) +
            scale_colour_manual(name = "Distributions",
                      labels = dists,
                      values = c("grey60", "black", "blue", "red")) +
            scale_linetype_manual(name = "Distributions",
                      labels = dists,
                      values = c("blank", "dashed", "dotted", "solid")) +
            scale_shape_manual(name = "Distributions",
                      labels = dists,
                      values = c(16, 26, 26, 26))

    leg <- g_legend(fig)
    dev.off()
    return(leg)
}


#######################################################################
#                    Plotting degree distributions                    #
#######################################################################


# plot degree distribution of the networks
plot_degree_distr <- function(ppi_name="string")
{
    # get the data
    data <- get_node_properties(ppi_name)

    # get the degrees
    degree <- as.integer(data$degree)

    degree_distr <- as.data.frame(prop.table(table(degree)))
    degree_distr$degree <- as.integer(levels(degree_distr$degree))


    # random graph (Erdos Renyi, or Gilbert)

    # get binomial random model
    binom_model <- fit_binom_degree_distr(degree)
    # cut of zero freq
    binom_model <- binom_model[which(binom_model$Freq > 1e-20),]

    # poisson model
    pois_model <- fit_poisson_degree_distr(degree)
    # cut of zero freq
    pois_model <- pois_model[which(pois_model$Freq > 1e-20),]

    # power law fit
    powerlaw_model <- fit_powerlaw_degree_distr(degree)
    # exponential distribution
    exp_model <- fit_exp_degree_distr(degree)


    fig <- ggplot(degree_distr, aes(x=degree, y=Freq)) + geom_point(colour="grey60")
    #fig <- fig + scale_x_log10() + scale_y_log10()
    fig <- fig + geom_line(data=binom_model, aes(x=degree, y=Freq), lty="dashed")
    #fig <- fig + geom_line(data=pois_model, aes(x=degree, y=Freq), lty="dotted")
    fig <- fig + geom_line(data=powerlaw_model, aes(x=degree, y=Freq), colour="red")
    fig <- fig + geom_line(data=exp_model, aes(x=degree, y=Freq), colour="blue", lty="dotted")
    fig <- fig + xlab("Degree") + ylab("Rel. Frequency")
    #fig <- fig + coord_cartesian(xlim=c(0,1000), ylim =c(min(degree_distr$Freq)/5, max(degree_distr$Freq) * 5))

    fig$degrees <- degree
    fig$degree_distr <- degree_distr

    return(fig)
}

# plot degree distribution of the networks
plot_degree_distr_no_fit <- function(ppi_name="string")
{
    # get the data
    data <- get_node_properties(ppi_name)

    # get the degrees
    degree <- as.integer(data$degree)

    degree_distr <- as.data.frame(prop.table(table(degree)))
    degree_distr$degree <- as.integer(levels(degree_distr$degree))


    fig <- ggplot(degree_distr, aes(x=degree, y=Freq)) + geom_point(colour="grey60")
    #fig <- fig + scale_x_log10() + scale_y_log10()
    fig <- fig + xlab("Degree") + ylab("Rel. Frequency")
    #fig <- fig + coord_cartesian(xlim=c(0,1000), ylim =c(min(degree_distr$Freq)/5, max(degree_distr$Freq) * 5))

    fig$degrees <- degree
    fig$degree_distr <- degree_distr

    return(fig)
}

plot_degree_distr_linear <- function(ppi_name="string")
{
    # get basic plot
    fig <- plot_degree_distr(ppi_name)
    degrees <- fig$degrees
    degree_distr <- fig$degree_distr

    # plots on linear scale (i.e. no log scale)
    adj_max <- quantile(degrees, 0.999) # throw away the top 1000th elements
    fig <- fig + coord_cartesian(xlim=c(0,adj_max), ylim =c(0, max(degree_distr$Freq)*1.05))

    return (fig)
}

plot_degree_distr_linear_no_fit <- function(ppi_name="string")
{
    # get basic plot
    fig <- plot_degree_distr_no_fit(ppi_name)
    degrees <- fig$degrees
    degree_distr <- fig$degree_distr

    # plots on linear scale (i.e. no log scale)
    adj_max <- quantile(degrees, 0.99) # throw away the top 1000th elements
    fig <- fig + coord_cartesian(xlim=c(0,adj_max), ylim =c(0, max(degree_distr$Freq)*1.05))

    return (fig)
}


plot_degree_distr_log <- function(ppi_name="string")
{
    # get basic plot
    fig <- plot_degree_distr(ppi_name)
    degrees <- fig$degrees
    degree_distr <- fig$degree_distr

    # plots on linear scale (i.e. no log scale)
    adj_max <- max(degrees) # throw away the top 1000th elements
    y_range <- c(min(degree_distr$Freq)/2, max(degree_distr$Freq) * 2)
    fig <- fig + scale_x_log10() + scale_y_log10()
    fig <- fig + coord_cartesian(xlim=c(1-0.,adj_max),ylim=y_range)

    return (fig)
}

plot_degree_distr_log_no_fit <- function(ppi_name="string")
{
    # get basic plot
    fig <- plot_degree_distr_no_fit(ppi_name)
    degrees <- fig$degrees
    degree_distr <- fig$degree_distr

    # plots on linear scale (i.e. no log scale)
    adj_max <- max(degrees) # throw away the top 1000th elements
    y_range <- c(min(degree_distr$Freq)/2, max(degree_distr$Freq) * 2)
    fig <- fig + scale_x_log10() + scale_y_log10()
    fig <- fig + coord_cartesian(xlim=c(1-0.,adj_max),ylim=y_range)

    return (fig)
}

#######################################################################
#                          string db example                          #
#######################################################################

plot_example_2 <- function(ppi_name="string")
{
    p <- ppi_name
    figs <- list()
    # plot linear plot (maybe without any distribution lines?)
    fig <- plot_degree_distr_linear(p)
    fig <- fig + labs(title="Degree Distribution of STRING")
    figs <- c(figs, list(fig))
    # TODO TODO

    # add legend??
    figs <- c(figs, list(degr_distr_legend()))

    all_figs <- do.call(grid.arrange, figs)
    return(all_figs)
}


plot_example_1 <- function(ppi_name="string")
{
    # simple linear plot of degree distribution
    p <- ppi_name
    figs <- list()
    fig <- plot_degree_distr_linear_no_fit(p)
    fig <- fig + labs(title="Degree Distribution of STRING")
    figs <- c(figs, list(fig))

    fig <- plot_degree_distr_log_no_fit(p)
    fig <- fig + labs(title="... with log-log scaling")
    figs <- c(figs, list(fig))

    all_figs <- do.call(grid.arrange, c(figs,list(nrow=1)))
    return(all_figs)
}

save_example_1 <- function()
{
    pdf("../figs/ppi_degr_distr.pdf", width=6.3, height=2.6)
    fig <- plot_example_1("string")
    print(fig)
    dev.off()
}


#######################################################################
#                          Helper functions                           #
#######################################################################


explode_degree_distr <- function(degree_distr)
{
    max_degr <- max(degree_distr$degree)
    expl_degree_distr <- data.frame(degree=1:max_degr, Freq=rep(0, max_degr))
    # fill the new distr with the old values where they fit
    for (i in 1:length(degree_distr$degree))
    {
        deg <- degree_distr$degree[i]
        frq <- degree_distr$Freq[i]
        expl_degree_distr$Freq[deg] = frq
    }
    return(expl_degree_distr)
}

expected_freq <- function(brks, cdf, cdf_params)
{
    est_freqs <- c()
    # for all intervals of the breaks
    for (i in 2:length(brks))
    {
        from <- brks[i-1]
        to <- brks[i]

        from_val <- do.call(cdf, c(from, cdf_params))
        to_val <- do.call(cdf, c(to, cdf_params))
        est_freqs <- c(est_freqs, to_val - from_val)
    }
    return(est_freqs)
}

test_fits <- function(ppi_name="string")
{
    # get the data
    data <- get_node_properties(ppi_name)

    # get the degrees
    degree <- as.integer(data$degree)

#    degree_distr <- as.data.frame(table(degree))
#    degree_distr$degree <- as.integer(levels(degree_distr$degree))
#
#    degree_distr <- explode_degree_distr(degree_distr)

    # TODO: unify the fits and tests into a single function
    # TODO: output table of fits per PPI (all are bad fits, but most are much
    #       much worse than powerlaw)
    # TODO: only fit for tail of degree distribution (i.e. ignoring low degree
    #       nodes in the evaluation of fit) [HOW? do i calc the freq for only the tail?]

    brks <- unique(quantile(degree, probs=seq(0,1, by=0.1)))
    c <- cut(degree, breaks=brks, include.lowest=TRUE)
    t <- as.data.frame(table(c))

    # random graph (Erdos Renyi, or Gilbert)

    # get binomial random model
    bin_expected_freq <- expected_freq(brks, pbinom, fit_binom(degree))
    #binom_model <- fit_binom_degree_distr(degree)
    # cut of zero freq

    fit <- chisq.test(t$Freq, p=bin_expected_freq, rescale.p=TRUE)
    print("Fit Binom")
    print(fit)


    # poisson model
    #pois_model <- fit_poisson_degree_distr(degree)
    bin_expected_freq <- expected_freq(brks, ppois, fit_poisson(degree))
    fit <- chisq.test(t$Freq, p=bin_expected_freq, rescale.p=TRUE)
    print("Fit Poisson")
    print(fit)
    # cut of zero freq

    # power law fit
    #powerlaw_model <- fit_powerlaw_degree_distr(degree)
    # TODO: this might fail because of xmin
    bin_expected_freq <- expected_freq(brks, ppowerlaw, fit_powerlaw(degree))
    fit <- chisq.test(t$Freq, p=bin_expected_freq, rescale.p=TRUE)
    #fit <- chisq.test(degree_distr$Freq, p=powerlaw_model$Freq, rescale.p=TRUE)
    print("Fit powerlaw")
    print(fit)

    # exponential distribution
    #exp_model <- fit_exp_degree_distr(degree)
    bin_expected_freq <- expected_freq(brks, pexp, fit_exp(degree))
    fit <- chisq.test(t$Freq, p=bin_expected_freq, rescale.p=TRUE)
    #fit <- chisq.test(degree_distr$Freq, p=exp_model$Freq, rescale.p=TRUE)
    print("Fit exp")
    print(fit)
}

#######################################################################
#                               all ppis distr                        #
#######################################################################


plot_all_ppis_degr_distr <- function()
{
    figs <- list()
    for (p in get_ppis())
    {
        fig <- plot_degree_distr_log(p)
        fig <- fig + labs(title=to_expr_name(p))
        figs <- c(figs, list(fig))
    }

    figs <- c(figs, list(degr_distr_legend()))

    all_figs <- do.call(grid.arrange, figs)
    return(all_figs)
}


save_all_ppi_degr_distr <- function()
{
    pdf("../figs/all_ppi_degr_distr.pdf", width=6, height=7)
    all_figs <- plot_all_ppis_degr_distr()
    print(all_figs)
    dev.off()
    return (figs)
}

# TODO: plot clustering coeff against p (N-1) [ see barabasi book ]
