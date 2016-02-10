# Copula transformation function (see Park and Gupta 2012), based on a variable's ECDF
	make_copula <- function(x) {
		return(ifelse(ecdf(x)(x)==1, qnorm(1-.0000001), qnorm(ecdf(x)(x))))
		}
