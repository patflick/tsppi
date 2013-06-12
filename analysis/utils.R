# Util functions
# 
# Author: Patrick Flick
###############################################################################


# a function that enables LHS list assignments
# 	source1: http://stackoverflow.com/questions/1826519/function-returning-more-than-one-value
# 	source2: https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
# enables assignment like:  list[a,,b] <- somefunction()
# which assigns a and b to the first and third element returned
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
	args <- as.list(match.call())
	args <- args[-c(1:2,length(args))]
	length(value) <- length(args)
	for(i in seq(along=args)) {
		a <- args[[i]]
		if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
	}
	x
}



