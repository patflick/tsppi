# plots for plotting the distrubtion of vertex properties of the two
# classes (for comparison). These are mostly general ways of comparing two
# distributions
# 
# Author: Patrick Flick
###############################################################################


plot_log_hist <- function(low_spec_properties, high_spec_properties, property_name)
{

    h_nbreaks = 7

    h_high <- hist(log10(high_spec_properties), breaks=h_nbreaks, plot=FALSE)
    h_low <- hist(log10(low_spec_properties), breaks=h_nbreaks, plot=FALSE)

    x_high <- 10**h_high$mids
    y_high <- h_high$counts/sum(h_high$counts)

    x_low <- 10**h_low$mids
    y_low <- h_low$counts/sum(h_low$counts)

    #plot(density(na.omit(high_spec_properties)), main="Density", col="red", log="xy")
    #lines(density(na.omit(low_spec_properties)), col="blue")
    plot(x_high, y_high, col="red", type="p", log="xy", xlim=c(1,max(x_high,x_low)), ylim=c(min(y_high, y_low),max(y_high,y_low)), pch=18, xlab=property_name, ylab="Fraction")
    lines(x_low, y_low, col="blue", type="p", pch=17)

    legend("topright",legend=c("Specific Proteins", "Promiscuous Proteins"), fill=c("red","blue"))

}


plot_log_hist_lm <- function(low_spec_properties, high_spec_properties, property_name)
{

    h_nbreaks = 7

    h_high <- hist(log10(high_spec_properties), breaks=h_nbreaks, plot=FALSE)
    h_low <- hist(log10(low_spec_properties), breaks=h_nbreaks, plot=FALSE)

    x_high <- 10**h_high$mids
    y_high <- h_high$counts/sum(h_high$counts)

    x_low <- 10**h_low$mids
    y_low <- h_low$counts/sum(h_low$counts)

    #plot(density(na.omit(high_spec_properties)), main="Density", col="red", log="xy")
    #lines(density(na.omit(low_spec_properties)), col="blue")
    plot(x_high, y_high, col="red", type="p", log="xy", xlim=c(1,max(x_high,x_low)), ylim=c(min(y_high, y_low),max(y_high,y_low)), pch=18, xlab=property_name, ylab="Fraction")
    abline(lm(log10(y_high)~log10(x_high)),col="red", lty=2)
    lines(x_low, y_low, col="blue", type="p", pch=17)
    abline(lm(log10(y_low)~log10(x_low)),col="blue",lty=2)


    legend("topright",legend=c("Specific Proteins", "Promiscuous Proteins"), fill=c("red","blue"))

}


plot_density <- function(low_spec_properties, high_spec_properties, property_name, log=FALSE)
{
    if (log)
    {
        plot(density(na.omit(high_spec_properties)), main="", col="red", log="xy", xlab=property_name, ylab="Density")
    }
    else
    {
        plot(density(na.omit(high_spec_properties)), main="", col="red", xlab=property_name, ylab="Density")
    }
    lines(density(na.omit(low_spec_properties)), col="blue")
    legend("topright",legend=c("Specific Proteins", "Promiscuous Proteins"), fill=c("red","blue"))
}





