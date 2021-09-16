# Credits: see https://cran.r-project.org/web/packages/Matrix/vignettes/Comparisons.pdf, for explanation
# on how to calculate OLS estimator in a fast way.

library(Matrix)
library(MASS)
require(compiler)
require(tidyr)

# transformation function for praise-winston auto-correlation correction
	praise_winsten <- function(x,rho) {
		res=double(length(x))
		xtrans = x[-1]
		xlag = x[-length(x)]
		
		res[1] = sqrt(1-rho^2)*x[1]
		res[-1] = xtrans-rho*xlag
		
		return(res)					
		}
	
	# compile function
	praise_winsten <- cmpfun(praise_winsten)

# iterative SUR
itersur <- function (X, Y, index, method = "FGLS", maxiter = 1000, reltol=10^-7, to.file=F) {
		
		if (!method%in%c('FGLS', 'FGLS-Praise-Winsten')) stop(paste('Invalid method selected: ',method))
		# verify correct data classes
		if (!class(X)=='matrix'|!class(Y)=='matrix') stop('X and Y need to be matrices')
		
		# verify order of index
		if (!all(order(index[,2], index[,1])==1:nrow(X))) stop('Data needs to be stacked by brands')
		colnames(index) <- c('date','brand')

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
				
		deltas <- NULL
		
		# iterate through SUR
		for (iter in 1:maxiter) {
			beta_old = beta_hat
			
			if (method=="FGLS-Praise-Winsten") {
				# Apply praise-winston correction for auto-correlation (see van Heerde, Leeflang, and Wittink, MktSci 2004, Greene 2003, section 20.9 (AR(1) disturbances)
				
				pred = X %*% beta_hat
				resid = Y - pred
				resid_by_brand = pivot_wider(data.frame(index, resid = matrix(resid)), names_from = "brand", values_from = "resid")
				rho_brand = apply(cbind(resid_by_brand[,-1]), 2, function(x) sum(x[-1]*x[-length(x)],na.rm=T)/sum(x^2,na.rm=T))
						
				yprime = matrix(unlist(mapply(praise_winsten, ysplit, as.list(rho_brand),SIMPLIFY=FALSE)),ncol=1)
				xprime = do.call('rbind', mapply(function(x,rho) apply(x, 2, praise_winsten, rho=rho), xsplit, as.list(rho_brand), SIMPLIFY = FALSE, USE.NAMES=FALSE))
				xprime=as(xprime, "dgeMatrix")
				}
			
			pred = xprime %*% beta_hat
			resid = yprime - pred
			
			resid_by_brand = pivot_wider(data.frame(index, resid = matrix(resid)), names_from = "brand", values_from = "resid")
			
			rhos = apply(cbind(resid_by_brand[,-1]), 2, function(x) sum(x[-1]*x[-length(x)],na.rm=T)/sum(x^2,na.rm=T))
			#print(rhos)
			sigma <- empty_sigma
			
			for (.i in 1:ncol(sigma)) {
				for (.j in 1:ncol(sigma)) {
					resids = cbind(resid_by_brand[, .i + 1], resid_by_brand[, .j + 1])
					compl.cases = complete.cases(resids)
					# if-statement below is true if no overlapping time periods
					if (length(which(compl.cases == TRUE)) <= 1) {
						sigma[.i, .j] <- 0
					}
					else {
						tmax = nrow(resids)
						resids = resids[complete.cases(resids), ]
						sigma[.i, .j] <- (1/tmax) * sum(resids[, 1] * resids[, 2])
					}
				}
			}
			
			sigma_inv = solve(sigma)
			#sigma = (1/tobs) * crossprod(as.matrix(resid_by_brand[,-1]))
			
			# old way to compute Kronecker, really slow
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
			
			inv_varcovar = crossprod(xprime, omega_inverse) %*% xprime
			varcovar = solve(inv_varcovar)
			beta_hat = varcovar %*% (crossprod(xprime, omega_inverse) %*% yprime)
			
			# check convergence, based on criterium in Greene (2002), p. 566
			delta = drop(t(beta_hat - beta_old) %*% inv_varcovar %*% (beta_hat - beta_old)) # the middle part belongs to the Hessian
			cat('Iteration ', iter, ' (Convergence Criteria: ', delta, ').\n')
			
			if (to.file==T) {
				outbetahat = data.frame(matrix(beta_hat,ncol=length(beta_hat)))
				
				colnames(outbetahat) <- colnames(xprime)
				outbetahat$iteration=iter
				outbetahat$delta=delta
				appendfile = ifelse(iter==1, F, T)
				write.table(outbetahat, 'iter_out.csv', append=appendfile, col.names = !appendfile, row.names=F)
				}
			
			deltas = c(deltas, delta)
			if (delta<reltol) break
		}
	
	res = new("itersur")
	
	pred = xprime %*% beta_hat
	resid = yprime - pred
	ses = sqrt(diag(varcovar))
	
	res@coefficients = data.frame(variable = colnames(X), coef = drop(as.matrix(beta_hat)), se = drop(as.matrix(ses)), ols = drop(as.matrix(beta_ols)), row.names=NULL)
	res@coefficients$z <- res@coefficients$coef/res@coefficients$se
	res@coefficients<-res@coefficients[, match(c('variable', 'coef', 'se', 'z'), colnames(res@coefficients))]
    
	k = ncol(varcovar)
    N = length(resid)
    res@bic = N * log(sum(resid^2)/N) + (k * log(N))
	res@llik = -.5 * N * log(sum(resid^2)/N)
	res@aic = 2*k - 2*res@llik
	
    res@predicted = as.numeric(pred)  # check with harald reg. computation of X, Y, etc.
    res@resid = as.numeric(resid)  # check with harald reg. computation of X, Y, etc.
	res@varcovar = as.matrix(varcovar)
    res@X <- as.matrix(xprime)
    res@y <- as.numeric(yprime)
    res@index <- index
	res@sigma <- as.matrix(sigma)
	res@iterations = iter
	res@delta = deltas
	res@method=method
	res@rho = rhos
	if(method=='FGLS-Praise-Winsten') res@rho_hat = rho_brand else res@rho_hat = rep(0, length(obsperbrand))
    return(res)
	
	# check with harald reg. computation of X, Y, etc.
	}



setClass("itersur",
	representation(X="matrix",
				   y="numeric",
				   index="data.frame",
				   sigma = "matrix",
				   iterations = "numeric",
				   delta = "numeric",
				   bic = "numeric",
				   aic = "numeric",
				   llik = "numeric",
				   predicted = "numeric",
				   resid = "numeric",
				   varcovar = "matrix",
				   coefficients = "data.frame",
				   method="character",
				   rho="numeric",
				   rho_hat="numeric"
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

			cat('\nEstimated auto-correlation of the residuals, if available:\n')
			print(object@rho_hat,digits=3)
			
			cat('\nRemaining auto-correlation of the residuals:\n')
			print(object@rho,digits=3)

			cat('\n')
			})
	
setMethod("coef", "itersur", function(object) {
			return(object@coefficients)
			
			})
# Wooldridge 2002; 8.36
# Rule of thumb: whenever residuals: use X; otherwise X hat.