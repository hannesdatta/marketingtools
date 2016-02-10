library(reshape2)
	sur <- function(X, Y, dates_brands) {
	
	# the columns of dates_brand needs to be: date, brand
	# the input need to be (stacked!): Ys, Xs.
	
	
	######################
	# ESTIMATE SUR MODEL #
	######################

		betas = solve(t(X)%*%X)%*%t(X)%*%Y
		pred = X%*%betas
		resid = Y-pred
		
		# Compute elements of covariance matrix s_hat (ij)
		resid_by_brand = dcast(data.frame(dates_brands, resid= matrix(resid)), date~brand,value.var='resid')
		resid_y_by_brand = data.table(data.frame(dates_brands, y=Y,resid= matrix(resid)))
		
		#r2ols = resid_y_by_brand[, list(r2=1-sum(resid^2)/sum((y-mean(y))^2)), by=c('brand')] #--> compute at a different point in time!
		
		 # r2 = 1-e'e/(y-meanc(y))'(y-meanc(y));
 
 
 
		# initialize var-covar matrix to zeros
		#if (0) {
		sigma <- matrix(double(length(unique(dates_brands$brand))^2),ncol=length(unique(dates_brands$brand)))

		# compute elements, see Bronnenberg, Mahajan, Vanhonacker 2000, p. 30, formula A2
		for (.i in 1:ncol(sigma)) {
			for (.j in 1:ncol(sigma)) {
				# i have to bring the residuals to one level
				resids = cbind(resid_by_brand[,.i+1],resid_by_brand[,.j+1])
				compl.cases = complete.cases(resids)
				if (length(which(compl.cases==T))<=1) {
					# no overlap; set sigma[.i, .j] to 0
					sigma[.i,.j] <- 0 } else {
					# overlap
					tmax = nrow(resids) # --> verify with Bart why it is tmax, and not tmin
					resids = resids[complete.cases(resids),]
					sigma[.i,.j]<- (1/tmax) * sum(resids[,1]*resids[,2])
					}
				}
			}
		#}
		
		#sigma = cov(resid_by_brand[,-1], use = 'pairwise.complete.obs')

		# in the derivations, I follow Greene (7th edition), p. 294, formula 10-7.
		max_t=length(unique(resid_by_brand$date))

		# obs per brand
		obsperbrand = apply(resid_by_brand[,-1], 2, function(x) length(which(!is.na(x))))
		
		sigma_inv = solve(sigma)
		
		#sigma_inv = backsolve(sigma, diag(ncol(sigma)))
		#sigma_inv[lower.tri(sigma_inv)] <- sigma_inv[upper.tri(sigma_inv)]
  
		# different observations per brand; need to compute Kronecker product manually
			inew = NULL
			for (.i in 1:ncol(sigma_inv)) {
				jnew = NULL
				for (.j in 1:ncol(sigma_inv)) {
				zeros = matrix(double(obsperbrand[.i]*obsperbrand[.j]),nrow=obsperbrand[.i], ncol=obsperbrand[.j])
				diag(zeros)<-sigma_inv[.i,.j]
				jnew = cbind(jnew, zeros)
				}
			inew = rbind(inew, jnew)
			}
		
		#omega = inew
		omega_inverse = inew # is sigma *kronecker* identiy matrix #solve(omega)

		#varcovar = backsolve(t(X) %*% omega_inverse %*% X, diag(ncol(X)))
		varcovar = solve(t(X) %*% omega_inverse %*% X)
		
		beta_hat = varcovar %*% (t(X) %*% omega_inverse %*% Y)
		
		data.frame(matrix(betas), matrix(beta_hat))
		
		
		# compute R2
		# residuals
		#resid = Y - X %*% beta_hat
		#resid_tmp = data.table(data.frame(dates_brands, y=Y,resid= matrix(resid)))
		#r2sur = resid_tmp[, list(r2=1-sum(resid^2)/sum((y-mean(y))^2)), by=c('brand')] 
		
		#
			
	res=NULL
	res$coefficients = data.table(variable=colnames(X), coef=drop(beta_hat), se = sqrt(diag(varcovar)))
	
	
	k=ncol(varcovar)
	N=length(resid)
	res$bic = log(sum(resid^2)/N) + (k*log(N))/N
	#res$r2 = r2
	#res$r2ols = r2ols
	#res$r2sur = r2sur
	res$predicted = X %*% beta_hat
	res$resid = Y - X %*% beta_hat
	res$X <- X
	res$Y <- Y
	res$dates_brands <- dates_brands
	#res$lags = lags
	return(res)
	
	#   e  = y - X*b;
   #ssr= e'e;
  # r2 = 1-e'e/(y-meanc(y))'(y-meanc(y));
  # AIC= ln(ssr/n) + (2*k)/n;
  # Schwarz = ln(ssr/n) + (k*ln(n))/n;

	}
	
