# Credits: see https://cran.r-project.org/web/packages/Matrix/vignettes/Comparisons.pdf, for explanation
# on how to calculate OLS estimator in a fast way.

library(Matrix)
library(MASS)

itersur <- function (X, Y, index, maxiter = 1000, method = "FGLS") {

		# verify correct data classes
		if (!class(X)=='matrix'|!class(Y)=='matrix') stop('X and Y need to be matrices')
		
		# verify order of index
		if (!all(order(index[,2], index[,1])==1:nrow(X))) stop('Data needs to be stacked by brands')
		
		X=as(X, "dgeMatrix")
		# initialize starting values for iterative SUR
			#beta_ols = solve(t(X) %*% X) %*% t(X) %*% Y
			beta_ols = solve(crossprod(X), crossprod(X, Y))
			beta_hat = beta_ols
		
		# get maximum observations by brand, used for Kronecker computations
			obsperbrand = table(index$brand)
			max_t = max(obsperbrand)
			takeouts <- c(sapply(obsperbrand, function(x) seq(from=1, to=max_t)%in%seq(from=1, to=x)))
			I=Diagonal(max_t)
		
		# initialize variance-covariance matrix
			empty_sigma <- sigma <- matrix(double(length(obsperbrand)^2), ncol = length(obsperbrand))
		
		# initialize x/y's-by-brand to ease computation for Praise-Winsten transformation
			xsplit = split(data.frame(as.matrix(X)), index[,2])
			ysplit = split(as.matrix(Y), index[,2])
			
			xprime = X
			yprime = Y
			
		# iterate through SUR
		for (i in 1:maxiter) {
			beta_old = beta_hat
			
			if (method=="FGLS-Praise-Winsten") {
				# Apply praise-winston correction for auto-correlation (see van Heerde, Leeflang, and Wittink, MktSci 2004, Greene 2003, section 20.9 (AR(1) disturbances)
				
				pred = X %*% beta_hat
				resid = Y - pred
				
				resid_by_brand = dcast(data.frame(index, resid = matrix(resid)), date ~ brand, value.var = "resid")
				rho_brand = apply(resid_by_brand[,-1], 2, function(x) sum(x[-1]*x[-length(x)])/sum(x^2))
						
				# apply transformation
				praise_winsten <- function(x,rho) {
					res=double(length(x))
					xtrans = x[-1]
					xlag = x[-length(x)]
					
					res[1] = sqrt(1-rho^2)*x[1]
					res[-1] = xtrans-rho*xlag
					
					return(res)					
					}
				
				yprime = matrix(mapply(praise_winsten, ysplit, as.list(rho_brand)),ncol=1)
				xprime = do.call('rbind', mapply(function(x,rho) apply(x, 2, praise_winsten, rho=rho), xsplit, as.list(rho_brand), SIMPLIFY =FALSE))
				xprime=as(xprime, "dgeMatrix")
				}
			
			pred = xprime %*% beta_hat
			resid = yprime - pred
			
			resid_by_brand = dcast(data.frame(index, resid = matrix(resid)), 
				date ~ brand, value.var = "resid")
			
			rhos = apply(resid_by_brand[,-1], 2, function(x) sum(x[-1]*x[-length(x)])/sum(x^2))
				
			sigma <- empty_sigma
			
			for (.i in 1:ncol(sigma)) {
				for (.j in 1:ncol(sigma)) {
					resids = cbind(resid_by_brand[, .i + 1], resid_by_brand[, 
						.j + 1])
					compl.cases = complete.cases(resids)
					if (length(which(compl.cases == TRUE)) <= 1) {
						sigma[.i, .j] <- 0
					}
					else {
						tmax = nrow(resids)
						resids = resids[complete.cases(resids), ]
						sigma[.i, .j] <- (1/tmax) * sum(resids[, 1] * 
						  resids[, 2])
					}
				}
			}
			
			sigma_inv = solve(sigma)
			#sigma = (1/tobs) * crossprod(as.matrix(resid_by_brand[,-1]))
			
			# old way
			if(0){
			inew = NULL
			for (.i in 1:ncol(sigma_inv)) {
				jnew = NULL
				for (.j in 1:ncol(sigma_inv)) {
					zeros = matrix(double(obsperbrand[.i] * obsperbrand[.j]), 
						nrow = obsperbrand[.i], ncol = obsperbrand[.j])
					diag(zeros) <- sigma_inv[.i, .j]
					jnew = cbind(jnew, zeros)
				}
				inew = rbind(inew, jnew)
			}
			omega_inverse = inew
			}
			
			# new way to compute Kronecker of unequal observations, much faster
			omega_inverse = kronecker(sigma_inv,I, make.dimnames = FALSE)[takeouts,takeouts]
			
			varcovar = solve(crossprod(xprime, omega_inverse) %*% xprime) #solve #ginv(as.matrix(crossprod(X, omega_inverse) %*% X)) #solve
	
			beta_hat = varcovar %*% (crossprod(xprime, omega_inverse) %*% yprime)
			
			# check convergence, based on criterium in Greene (2002), p. 566
			delta = drop(t(beta_hat - beta_old) %*% (t(xprime)%*%omega_inverse%*%xprime) %*% (beta_hat - beta_old)) # the middle part belongs to the Hessian
			if (delta<10^-7) break
		}
		
	res = new("itersur")
	
	varcovar = solve(t(xprime)%*%omega_inverse%*%xprime)
	ses = sqrt(diag(varcovar))

	res@coefficients = data.frame(variable = colnames(X), coef = drop(as.matrix(beta_hat)), se = drop(as.matrix(ses)), ols = drop(as.matrix(beta_ols)), row.names=NULL)
	res@coefficients$z <- res@coefficients$coef/res@coefficients$se
	res@coefficients<-res@coefficients[, match(c('variable', 'coef', 'se', 'z'), colnames(res@coefficients))]
    k = ncol(varcovar)
    N = length(resid)
    res@bic = log(sum(resid^2)/N) + (k * log(N))/N
    res@predicted = as.numeric(xprime %*% beta_hat)  # check with harald reg. computatoin of X, Y, etc.
    res@resid = as.numeric(yprime - xprime %*% beta_hat)  # check with harald reg. computatoin of X, Y, etc.
	res@varcovar = as.matrix(varcovar)
    res@X <- as.matrix(xprime)
    res@y <- as.numeric(yprime)
    res@index <- index
	res@sigma <- as.matrix(sigma)
	res@iterations = i
	res@method=method
	res@rho = rhos
    return(res)
	
	# check with harald reg. computatoin of X, Y, etc.
	}



setClass("itersur",
	representation(X="matrix",
				   y="numeric",
				   index="data.frame",
				   sigma = "matrix",
				   iterations = "numeric",
				   bic = "numeric",
				   predicted = "numeric",
				   resid = "numeric",
				   varcovar = "matrix",
				   coefficients = "data.frame",
				   method="character",
				   rho="numeric"
				   ))

setMethod("show", "itersur", function(object) {
			options(scipen=999)
			cat('itersur() Estimation (FGLS) \n')
			cat('=============================\n\n')
			
			cat('Estimation method              : ', object@method,'\n')
			cat('Iterations                     : ', object@iterations,'\n')
			cat('BIC:                           : ', object@bic, '\n')
			
			cat('\nParameter Estimates:\n')
			print(object@coefficients,digits=3)
			cat('\n')
			
			cat('\nCross-sectional variance-covariance matrix:\n')
			print(object@sigma,digits=3)
			cat('\n')

			cat('\nAuto-correlation of the residuals:\n')
			print(object@rho,digits=3)
			cat('\n')
			})
	
setMethod("coef", "itersur", function(object) {
			return(object@coefficients)
			
			})
