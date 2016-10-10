require(data.table)
require(reshape2)	
require(Formula)	


attraction_data <- function(formula, data = NULL, heterogenous=NULL, index = NULL, model = 'MCI', benchmark = NULL, ...) { #subset
	# converts a data set to its corresponding attraction-based data (using the base-brand approach), see Fok 2001.
	# index: first column is individual, second is time

	# Retrieve all variables
		xhom = cbind(stats::model.frame(formula=formula, data=data,na.action=NULL)[,-1])
		y = stats::model.frame(formula=formula, data=data,na.action=NULL)[,1]
		xhet = cbind(stats::model.frame(formula=heterogenous, data=data,na.action=NULL))
		cnames_xhom = colnames(xhom)[!colnames(xhom)%in%colnames(xhet)]
		
		xhom <- cbind(xhom[,!colnames(xhom)%in%colnames(xhet)])
		colnames(xhom) <- cnames_xhom
		index = stats::model.frame(formula=index, data=data,na.action=NULL)
		
		if (!ncol(index)==2) stop("Please provide an index: brands (column 1), time (column 2)")
	
	n_individ = length(unique(index[,1]))
	
	# create conversion matrix for base-brand approach
	Hbb = diag(n_individ-1)
	Hbb = cbind(Hbb, -1)

	# Transformed X and Y matrices
	dtmelt <- reshape2::melt(data.frame(y=y, xhom, xhet, individ=index[,1], period=index[,2]), id.vars=c('individ', 'period'))
	
	# select benchmark brand: the one with most available observations (mean of all variables)
	tmp <- table(dtmelt[which(dtmelt$variable=='y' & !is.na(dtmelt$value)), c('individ')])
	if (is.null(benchmark)) {
		bindivid = rev(names(tmp)[which(tmp==max(tmp))])[1] } else {
		bindivid = benchmark }
		
	aindivid = names(tmp)[!names(tmp)==bindivid]
	
	# stacked data set
	dtcast = data.table::rbindlist(lapply(aindivid, function(x) {
		tmp=reshape2::dcast(dtmelt[dtmelt$individ%in%c(x, bindivid),], period ~ individ + variable)
		# kick out variables which are unavailable throughout the whole period
		tmp=tmp[,!colSums(is.na(tmp))==nrow(tmp)]
		tmp[complete.cases(tmp),]
		}), fill=TRUE)

	iindex = aindivid[apply(dtcast[, paste0(aindivid,'_y'),with=F],1, function(x) which(!is.na(x)))]
	
	# apply transformations
		if (model=='MCI') tfkt <- function(x) log(x)
		if (model=='MNL') tfkt <- function(x) x
	
	# market share transformations!
	
	# 
	
	# heterogenous transformation:
	# find variables which are heterogenous
	dttrans <- dtcast
	
	# transformation for heterogenous parameters
	for (.var in grep(paste(paste0('[_]', colnames(xhet)),collapse='|'), colnames(dtcast),value=T)) {
		if (grepl(paste0('^',bindivid),.var)) {
			# variable pertaining to base brand
			dttrans[, .var:= -tfkt(get(.var)),with=F]
			} else {
			# variable pertaining to other brands
			dttrans[, .var:= tfkt(get(.var)),with=F]
			}
		}
		
	# transformation of 'homogenous' coefficients
	for (.var in colnames(xhom)) {
		# identify target cols
		tmp1=rowSums(dttrans[, mget(paste0(aindivid, '_', .var))], na.rm=T)
		tmpbase=dttrans[, mget(paste0(bindivid, '_', .var))]
		
		dttrans[, paste0('hom_', .var) := tfkt(tmp1) - tfkt(tmpbase),with=F]
		dttrans[, c(paste0(bindivid, '_', .var),paste0(aindivid, '_', .var)):=NULL, with=F]
		}

	# transformation of 'y' coefficients
	for (.var in c('y')) {
		# identify target cols
		tmp1=rowSums(dttrans[, mget(paste0(aindivid, '_', .var))], na.rm=T)
		tmpbase=dttrans[, mget(paste0(bindivid, '_', .var))]
		
		dttrans[, paste0('hom_', .var) := log(tmp1/tmpbase),with=F]
		dttrans[, c(paste0(bindivid, '_', .var),paste0(aindivid, '_', .var)):=NULL, with=F]
		}

	dummies = as.matrix(model.matrix( ~ as.factor(iindex) - 1))
	colnames(dummies)<-paste0(aindivid,'_dum')
	
	dttrans[is.na(dttrans)]<-0
	
	out = new("attr.data")
	# remove homogenous (pre-transformation columns) from dttrans
	xmatrix = as.matrix(cbind(dttrans[, !grep('hom_y|period', colnames(dttrans),value=T),with=F], dummies))
	#rownames(xmatrix) <- rep(aindivid, 
	out@X <- xmatrix	
	out@y <- as.numeric(dttrans$hom_y)
	out@period <- as.numeric(dttrans$period)
	out@individ <- iindex
	out@input <- list(formula=formula, data=data, index=index)
	out@model <- model
	out@benchmark <- bindivid
	
	return(out)
	}



setClass("attr.data",
	representation(X="matrix",
				   y="numeric",
				   individ = "character",
				   period = "numeric",
				   model = "character",
				   benchmark = "character",
				   input = "list"),
	
	validity = function(object) {
		retval=NULL
		
		if (!length(unique(colnames(object@X)))==ncol(object@X)) retval <- c(retval, 'Column names of X not properly specified')
		if (ncol(object@X)==0|nrow(object@X)==0) retval <- c(retval, 'X has 0-dimension')
		
		if (is.null(retval)) return(TRUE) else return(retval)
		})



setMethod("show", "attr.data", function(object) {
			cat('Attraction model: Data set\n')
			
			
			cat('    Model type                                    : ', object@model,'\n')
			cat('    Benchmark brand name                          : ', object@benchmark,'\n\n\n')
			
			cat('    Number of total cross-sectional units         : ', length(unique(object@input[["index"]][,1])), '\n')
			cat('    Number of total periods                       : ', length(unique(object@input[["index"]][,2])), '\n\n')
			
			cat('    Number of variables in transformed system     : ', ncol(object@X), '\n')
			cat('    Observations per cross-sectional unit         : ', table(object@individ), '\n')
			cat('    Column names of X                             : ',paste(colnames(object@X), collapse= ', '),'\n')
			cat('\n')
			
			})




