attraction

# Simulate some data

makedata <- function() {
		
	# Simulate sample data

	n_brands = 3
	set.seed(1984)

	tint = c(-1,-2,0) #c(runif(n_brands-1),0) # true intercepts, of which the first one should be zero always

	tbeta = matrix(c(.5, 1, # true betas (note: heterogenous coefficients)
					 0, 1,
					-1,-1), ncol=2, byrow=T)

	tattr = c(.5,1,-.5)

	tobs = 1000 # number of observations

	xsim = matrix(runif(tobs*(length(tbeta)+length(tattr)*n_brands)),ncol=ncol(tbeta)+length(tattr))
	pars_heterog = 1:ncol(tbeta)
	pars_homog = (ncol(tbeta)+1):ncol(xsim)

	# set benchmark brand to last brand in this simulation.
	benchm <- n_brands

	sigma_tilde <- matrix(c(1,.25, .25, .5), ncol=2)

	if (!length(unique(c(n_brands, length(tint), ncol(sigma_tilde)+1, nrow(sigma_tilde)+1)))==1) stop('Check true parameters; some problem with dimensions.')
	if (!nrow(tbeta)==n_brands) stop('Heterogenous coefficients not specified properly')

	set.seed(1984)

	y <- double(tobs*n_brands)

	for (p in 1:tobs) { # time periods
		nu = mvrnorm(n=1, mu=rep(0, ncol(sigma_tilde)), Sigma = sigma_tilde)
		
		index = seq(from=p*n_brands-(n_brands-1), length.out=n_brands)
		
		tbetaattr <- cbind(tbeta, matrix(rep(tattr, n_brands),ncol=length(tattr),byrow=T))
		
		xbeta = xsim[index,] ^ rbind(tbetaattr[-benchm,], -tbetaattr[benchm,])
		
		# compute m_it
		m <- double(n_brands)
		m[benchm] <- 1
		m[-benchm] <- exp(tint[-benchm]+nu) * 
					  apply(xbeta[-benchm,], 1, "prod") * 
					  rep(prod(xbeta[benchm,]), n_brands-1)
		
		M=m/sum(m)
		
		y[index] <- M
		}

	index = data.frame(var=as.factor(rep(1:(n_brands), tobs)), t=rep(1:tobs,each=n_brands))

	X<<-xsim
	y<<-y
	indexv<<-index
	}
	
	
# Introduction to making classes in R

setClass("attr.raw",
	representation(X="matrix",
				   y="numeric",
				   individ = "character",
				   period = "numeric"),
	validity = function(object) {
		retval=NULL
		
		if (!length(unique(colnames(object@X)))==ncol(object@X)) retval <- c(retval, 'Column names of X not properly specified')
		if (ncol(object@X)==0|nrow(object@X)==0) retval <- c(retval, 'X has 0-dimension')
		
		if (is.null(retval)) return(TRUE) else return(retval)
		})
				   
new("attr.raw")

setClass("attr.bb", representation(X="matrix",
								   y="numeric",
								   individ = "numeric",
								   period = "numeric"),
					contains = "attr.raw")


					
dt = new("attr.raw")
validObject(dt)

makedata()

colnames(X) <- paste0('var_', 1:ncol(X))
dt@X <- as.matrix(X)
dt@y <- as.numeric(y)
dt@individ <- as.character(indexv$var)
dt@period <- indexv$t
validObject(dt)

setMethod("show", "attr.raw", function(object) {
			cat('Attraction model: Summary of raw data\n')
			cat('    Number of total cross-sectional units         : ', length(unique(object@individ)), '\n')
			cat('    Number of total periods                       : ', length(unique(object@period)), '\n')
			cat('    Number of variables                           : ', ncol(object@X), '\n')
			cat('    Observations per cross-sectional unit         : ', table(object@individ), '\n')
			cat('    Column names of X                             : ',paste(colnames(object@X), collapse= ', '),'\n')
			cat('\n')
		
			# check number of overlapping observations
			
			
			})
			
show(dt)

setGeneric("convertbb", function(object, ...) {})

require(data.table)

setMethod("convertbb", "attr.raw", function(object, model = 'MCI', heterogenous = rep(0, ncol(object@X))) {
	#print(heterogenous)
	#print(object)
	n_individ = length(unique(object@individ))
	
	# Fixes base-brand to be the last brand (for prototyping)
	
	# create conversion matrix for base-brand approach
	Hbb = diag(n_individ-1)
	Hbb = cbind(Hbb, -1)

	# Transformed X and Y matrices
	#x_hom <- object@X[,heterogenous==0]
	#x_het <- object@X[,heterogenous==1]
	
	dtmelt <- melt(data.frame(object@X, y=object@y, individ=object@individ, period=object@period), id.vars=c('individ', 'period'))
	
	# select benchmark brand: the one with most available observations (at a tie, take last one)
	tmp <- table(dtmelt[which(dtmelt$variable=='y' & !is.na(dtmelt$value)), c('individ')])
	bindivid = rev(names(tmp)[which(tmp==min(tmp))])[1]
	aindivid = names(tmp)[!names(tmp)==bindivid]
	
	# stacked data set
	dtcast = rbindlist(lapply(aindivid, function(x) {
		tmp=dcast(dtmelt[dtmelt$individ%in%c(x, bindivid),], period ~ individ + variable)
		tmp[complete.cases(tmp),]
		}), fill=TRUE)
	
	
	iindex = aindivid[apply(dtcast[, mget(paste0(aindivid,'_y'))],1, function(x) which(!is.na(x)))]
	
	#dtcast[, individ:=iindex]
	
	# apply transformations
	#if (model=='MCI')
	
	# heterogenous transformation:
	# find variables which are heterogenous
	x_hom <- colnames(object@X)[which(heterogenous==0)]
	x_het <- colnames(object@X)[which(heterogenous==1)]
	
	x_hom2 <- grep(paste(x_hom,collapse='|'), colnames(dtcast),value=T)
	x_het2 <- grep(paste(x_het,collapse='|'), colnames(dtcast),value=T)
	y_hom2 <- grep(paste('y',collapse='|'), colnames(dtcast),value=T)
	
	dttrans <- dtcast
	
	# transformation of 'heterogenous' coefficients
	for (.var in x_het2) {
		if (grepl(paste0('^',bindivid),.var)) {
			# variable pertaining to base brand
			dttrans[, .var:= -log(get(.var)),with=F]
			} else {
			# variable pertaining to other brands
			dttrans[, .var:= log(get(.var)),with=F]
			}
		}
		
	# transformation of 'homogenous' coefficients
	for (.var in c('y', x_hom)) {
		# identify target cols
		tmp1=rowSums(dttrans[, mget(paste0(aindivid, '_', .var))], na.rm=T)
		tmpbase=dttrans[, mget(paste0(bindivid, '_', .var))]
		
		dttrans[, paste0('hom_', .var) := log(tmp1/tmpbase),with=F]
		dttrans[, c(paste0(bindivid, '_', .var),paste0(aindivid, '_', .var)):=NULL, with=F]
		}
	
	dummies = as.matrix(model.matrix( ~ as.factor(iindex) - 1))
	colnames(dummies)<-paste0('dum_',aindivid)

	
	# remove homogenous (pre-transformation columns) from dttrans
	dttrans[is.na(dttrans)]<-0
	object@X <- as.matrix(cbind(dttrans[, !grep('hom_y|period', colnames(dttrans),value=T),with=F], dummies))
	object@y <- as.numeric(dttrans$hom_y)
	object@period <- as.numeric(dttrans$period)
	object@individ <- iindex
	return(object)
	})
	
dt2 <- convertbb(dt, heterogenous = c(1,1,0,0,0))


