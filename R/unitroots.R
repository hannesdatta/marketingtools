#############################################
### Functions to conduct unit root tests  ###
#############################################

library(fUnitRoots) # used to obtain McKinnon's (1996) p-values for Unit Root Tests

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
		
		list(t=tstat, p = pval,fstatistic = sum.m$fstatistic, SBC = bic, SSR=sum(resid(m)^2),n=n,k=k,m1=m)
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
	
#adf_enders_single(x$used_value, maxlag= max.lag, pval = pval, season=NULL)
#		X=x$used_value
		
# Conducts ADF test using Elder and Kenny's procedure on X (which can be either differenced, or non-differenced)
	adf_enders_single <- function(X, maxlag, season, pval) {
		eq1=adf.test(X, maxlag=maxlag, trend=T,season=season)
		testpval = eq1$p
		
		# Procedure described in Elder and Kennedy: Testing for Unit Roots: What Students Should be Taught?
		if (testpval<=pval) {  
			# --> rejection of test: no unit root; but: test for trend
			
			# next: test for significance of trend
			if ((summary(eq1$m1)$coefficients["lintrend",4])<=pval) {
				# trend significant --> NO UR, and a deterministic trend
				res = list(trend=1, ur = 0, order = 0,n=eq1$n, t=eq1$t,p=testpval, lags = eq1$lags)
				} else {
				# trend insignificant --> NO UR, no deterministic trend
				res = list(trend=0, ur = 0, order = 0,n=eq1$n, t=eq1$t,  p=testpval, lags = eq1$lags)
				}
			} else {
			
			# --> not reject: unit root, but NO deterministic trend (check for higher orders!)
			res = list(trend=0, ur = 1, order = 1, n=eq1$n, t=eq1$t, p=testpval, lags = eq1$lags)
			}
		return(res)
		}

	#adf_regular <- function(X, maxlag, season=NULL) {
	#	eq1=adf.test(X, maxlag=maxlag, trend=F,season=season)
	#	testpval = punitroot(eq1$t, N=eq1$n, trend = 'c', statistic='t')
	#	return(testpval)
	#	}
	#adf_regular<-cmpfun(adf_regular)
	
	
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

	
	
# combine p values, see Verbeke 2004, p. 372
	#adf_by_brand[, ':=' (N_brands=.N, P=-2*sum(log(p))),by=c('category','country', 'variable')][order(category, country,variable)]
	#adf_by_brand[, combined_p_val := pchisq(unique(P), df=2*unique(N_brands), lower.tail=FALSE), by=c('category','country', 'variable')]
	#setorderv(adf_by_brand, c('category', 'country','variable'))
