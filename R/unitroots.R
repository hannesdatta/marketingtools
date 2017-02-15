#############################################
### Functions to conduct unit root tests  ###
#############################################

library(fUnitRoots) # used to obtain McKinnon's (1996) p-values for Unit Root Tests
library(data.table) # to use the shift function

# Function to difference a series
	makediff = function(x, order=1) {
		res=x
		for (i in seq(length.out=order)) {
			res=res-makelag(res)
			}
		return(res)
		}
	
	
# Function to lag a series		
	makelag = function(x, lags = 1) {
		c(rep(NA,lags),x)[1:length(x)]
		}

		
# Function to conduct UR tests concurrently for up to maxlag lags using the ADF test equation with a constant and (optional) trend

adf.test <- function(X, maxlag = 12, trend = TRUE, season = NULL) {
	#maxlag = 12 # number of maximum lags when determining "common" size of test (to make sure BICs are computed even for 1 lag for data in which lags 12 are used)
	#lags = 1:12 # lags to be tested
	#trend = T
	##season is not activated
	if (!is.null(season)) stop('Seasonality not supported yet')
	
	# Create variables for ADF test equation:
	# delta_Y_t = alpha + beta*trend + gamma*y_t-1 + d_1 * delta_y_t-1 + d_2 * delta_y_t-2, etc.
	# test on coefficient of gamma
	
	z_diff = c(NA,diff(X))
	lag_matrix = matrix(unlist(shift(z_diff,1:maxlag)),ncol=maxlag)
	z_lag_1 = shift(X, 1)

	# Determine first and last observation in data
	non_na = which(!is.na(X))
	first_nonna = non_na[1]
	last_nonna = rev(non_na)[1]

	actual_maxlag = pmin(maxlag, last_nonna - first_nonna - 1)
	
	# Determine observations to estimate models on (i.e., all lags definining max lag are taken out)
	subs = seq(from = first_nonna + actual_maxlag + 1, to = last_nonna)
		
	# Add linear trend, if necessary
	lintrend = NULL
	if (trend==T) lintrend = subs-min(subs)+1

	# Run up to 0 to maxlag regressions, and save bic to determine selected lag length
	y = z_diff[subs]
	n = length(y)
	
	bics = sapply(0:actual_maxlag, function(lags) {
		.seq=seq(to=ncol(lag_matrix)-(actual_maxlag-lags), length.out=lags)
		df = as.matrix(cbind(lag1 = z_lag_1[subs], intercept=1, lintrend = lintrend, lag_matrix[subs, .seq]))
		XprimeXinv = solve(t(df)%*%df)
		beta=XprimeXinv%*%(t(df)%*%y)
		resid = y-df%*%beta
		k = ncol(df)
		bic = n * log(sum(resid^2)/n) + (k * log(n))
		})

	# Define selected lag
	sellag = which(bics==min(bics))-1

	# Re-estimate final model on all possible observations
	subs = seq(from = first_nonna + sellag + 1, to = last_nonna)
	lintrend = NULL
	if (trend==T) lintrend = subs-min(subs)+1

	y = z_diff[subs]
	n = length(y)

	.seq=seq(to=ncol(lag_matrix)-(maxlag-sellag), length.out=sellag)
	df = as.matrix(cbind(lag1 = z_lag_1[subs], intercept=1, lintrend = lintrend, lag_matrix[subs, .seq]))
			
	XprimeXinv = solve(t(df)%*%df)
	beta=XprimeXinv%*%(t(df)%*%y)
	rownames(beta)[which(rownames(beta)=='')] <- paste0('lag_delta_t-', seq(1, length.out = sellag))
	
	resid = y-df%*%beta
	k = ncol(df)
	varcovar = 1/(n-k) * drop(t(resid)%*%resid) * XprimeXinv
	se = sqrt(diag(varcovar))
	
	tstat=beta[1]/se[1]
	if(trend==T) trendspec = 'ct' else trendspec ='c'
	pval = fUnitRoots::punitroot(tstat, N=n, trend = trendspec, statistic='t')
	
	names(bics) <- paste0('lag_length_', 0:(length(bics)-1))
	
	coefs = data.frame(beta=beta, se=se)
	coefs$t = coefs$beta/coefs$se
	coefs$p = (1-pnorm(abs(coefs$t)))*2
	coefs$p = (1-pt(abs(coefs$t), n-k-1))*2

	trend_p = NA
	if (trend==T) trend_p = coefs[3,4]
	
	list(t=tstat, p=pval, trend_p=trend_p, coefs = cbind(beta, se), bics = bics, lags = sellag, n = n, k = k, trend = trend)
	}

if(0){
# Function to conduct unit root tests with fixed number of lags (on a sample restricted on observations starting at maxlag, to keep up the 
	adf.equation <- function(X, lags, maxlag, trend, season) { # restrict estimation to elements 'as if' there were ten lags (used for comparison)
		z.diff = makediff(X)
		lag_matrix = cbind(sapply(0:lags, function(l) makediff(makelag(X,l)))[,-1])
		lag_matrix_virt = cbind(sapply(0:maxlag, function(l) makediff(makelag(X,l)))[,-1])
		z.lag.1 = makelag(X,1)
		if (is.null(maxlag)) maxlag=lags
		
		df = data.frame(z.diff, lag_matrix, z.lag.1)
		
		subs = complete.cases(df, lag_matrix_virt) # determine on which elements the ADF model will be estimated (accoutning for different no. of lagged terms)
		#table(subs)
		
		lintrend = rep(NA, length(X))
		lintrend[subs] <- seq(from=1,length.out=length(which(subs==T)))
		
		form = z.diff~1+ z.lag.1
		if (trend==T) form=update(form, ~. + lintrend)
		if (length(season)>0 & length(unique(season)>1)) form=update(form, ~. + as.factor(season))
		if (ncol(lag_matrix)>0) form=update(form, ~ . + lag_matrix)

		m<-lm(form, subset=subs)
		sum.m = summary(m)
		#summary(m)
		
		n=length(which(!is.na(z.diff[subs])))
		k=length(m$coefficients)
		bic = (as.numeric(logLik(m))*-2)/n+(k*log(n))/n
		
		if(trend==T) trendspec = 'ct' else trendspec ='c'
		tstat=sum.m$coefficients["z.lag.1",3]
		pval = fUnitRoots::punitroot(tstat, N=n, trend = trendspec, statistic='t')
		
		list(t=tstat, p = pval, fstatistic = sum.m$fstatistic, SBC = bic, SSR=sum(resid(m)^2),n=n,k=k,m1=m)
		}

	
# conducts ADF test (as specified in adf.equation) maxlag times, and selects best-fitting stats
	adf.test <- function(X, maxlag, trend, season) {
		
		# compute best-fitting BIC per model
		bics=unlist(sapply(0:maxlag, function(i) adf.equation(X, lags=i, maxlag=maxlag,trend=trend,season=season)$SBC))
		# select best-fitting model
		sellag = which(bics==min(bics))-1
		
		# reestimate with selected lags
		out=adf.equation(X, lags=sellag, maxlag=sellag, trend=trend,season=season)
		out$lags<-as.integer(sellag)
		return(out)
		}
	
	
}
	
#adf_enders_single(x$used_value, maxlag= max.lag, pval = pval, season=NULL)
#		X=x$used_value
		
# Conducts ADF test using Elder and Kenny's procedure on X (which can be either differenced, or non-differenced)
	adf_enders_single <- function(X, maxlag, season, pval) {
		eq1=adf.test(X, maxlag=maxlag, trend=T, season=season)
		
		names(eq1$t) <- NULL
		names(eq1$p) <- NULL
			
		# Procedure described in Elder and Kennedy: Testing for Unit Roots: What Students Should be Taught?
		if (eq1$p<=pval) {
			# --> rejection of test: no unit root; but: test for trend
			
			
			# next: test for significance of trend
			if (eq1$trend_p<=pval) {
				# trend significant --> NO UR, and a deterministic trend
				res = list(trend=1, ur = 0, order = 0, n=eq1$n, t=eq1$t, p=eq1$p, lags = eq1$lags)
				} else {
				# trend insignificant --> NO UR, no deterministic trend
				res = list(trend=0, ur = 0, order = 0,n=eq1$n, t=eq1$t,  p=eq1$p, lags = eq1$lags)
				}
			} else {
			# --> not reject: unit root, but NO deterministic trend (check for higher orders!)
			res = list(trend=0, ur = 1, order = 1, n=eq1$n, t=eq1$t, p=eq1$p, lags = eq1$lags)
			}
		return(res)
		}

###########
## MAIN  ##
###########

# wrapping function to execute Elders and Kenny's test for unit roots
adf_enders <- function(X, ...) {
	# check number of unique values
	if ((length(unique(X))<2)|all(is.na(X))) return(NULL)
	
	difforder = 0:2
	res=do.call('rbind',lapply(difforder, function(.lag) adf_enders_single(makediff(X,.lag), ...)))
	
	out=res[1,]
	out$order=0+out$ur
	
	if ((res[1,]$ur==1 & res[2,]$ur==1)) {
		out=res[2,]
		out$order = 2
		}

	if ((res[1,]$ur==1 & res[2,]$ur==1 & res[3,]$ur==1)) {
		out=res[3,]
		out$order = 3
		}
	return(do.call('c',out))
		
	}

#	adf.test(X)


# Test:
# adf_enders(X, maxlag = 12, season = NULL, pval = .05)


# combine p values, see Verbeke 2004, p. 372
	#adf_by_brand[, ':=' (N_brands=.N, P=-2*sum(log(p))),by=c('category','country', 'variable')][order(category, country,variable)]
	#adf_by_brand[, combined_p_val := pchisq(unique(P), df=2*unique(N_brands), lower.tail=FALSE), by=c('category','country', 'variable')]
	#setorderv(adf_by_brand, c('category', 'country','variable'))
