source('attraction.R')
source('itersur.R')
source('attraction_data.R')


# Simulate data
	rawdata <- do.call('cbind', attraction_simulate_data())

# Define model
	head(rawdata)
	formula = y ~ X.var_3 + X.var_4 + X.var_5 | X.var_1 + X.var_2 | index.var + index.t
	subset = NULL
	benchmark=NULL

# Transform to base-brand system of equations representation
	dtbb <- attraction_data(formula, data=rawdata, model="MCI")
	validObject(dt)
	show(dt)

# Estimate the system
	m <- itersur(X=dtbb@X,Y=as.matrix(dtbb@y), index=data.frame(date=dtbb@period,brand=dtbb@individ))
	m
	coef(m)


#	m$coefficients$variable <- c(paste0(rep(paste0('var', 1:ncol(xsim[, pars_heterog]), '_'), each=n_brands), rep(paste0('b', 1:n_brands), ncol(xsim[, pars_heterog]))),
#								 paste0(rep(paste0('attr', 1:ncol(xsim[, pars_homog]), '_'), each=1)), 
#								 paste0('dum', 1:(n_brands-1)))

	#m$coefficients$true=c(matrix(tbeta, nrow=length(unlist(tbeta)), byrow=T), tattr, tint[-n_brands])

	print(m$sigma)
	m$coefficients

# HOW ARE LAGGED DVs INCLUDED?

# Include them as heterogenous effects


	
# Next steps: calculate elasticities
	m <- itersur(X=dat@X,Y=as.matrix(dat@y), dates_brands=data.frame(date=dat@period,brand=dat@individ))
	
	m$coefficients$variable <- c(paste0(rep(paste0('var', 1:ncol(xsim[, pars_heterog]), '_'), each=n_brands), rep(paste0('b', 1:n_brands), ncol(xsim[, pars_heterog]))),
								 paste0(rep(paste0('attr', 1:ncol(xsim[, pars_homog]), '_'), each=1)), 
								 paste0('dum', 1:(n_brands-1)))

	ols=solve(t(X)%*%X)%*%t(X)%*%Y
	#m$coefficients$true=c(matrix(tbeta, nrow=length(unlist(tbeta)), byrow=T), tattr, tint[-n_brands])
	m$coefficients$ols=ols

	print(m$sigma)
	m$coefficients


# Simulation



# Implementing the model formula object

source('attraction.R')
source('itersur.R')
source('attraction_data.R')


# The first describes possible Xs and possible Ys; it is verifying the availability of all Xs and all Ys, and - if necessary - kicks out observations that are not complete.
# There is no need to specify whether a given model is homogenous or heterogenous; as transformations will only be applied later.

# Data prep 1: verify structure of data set
#  o List possible base brands
#  o Verify data availability
#  o 

# Data prep 2: convert to BB approach
#  o Explicitly specify how many lagged values are to be included



# Formula: dependent variable (e.g., a market share) ~ homogenous (e.g., 1, an intercept), and | heterogenous coefficients (dummy variables cannot be done here); the homogenous OR heterogenous component should be able to contain lagged depedent variables.



require(data.table)
require(reshape2)



	
	# compute market shares
	mf[, ms := y/sum(y), by = c('index_time')]
	setorder(mf, index_individ, index_time)
	
	
	# Populate dt
	dt@X <- as.matrix(rawdata$X)
	dt@y <- as.numeric(rawdata$y)
	dt@individ <- as.character(rawdata$index$var)
	dt@period <- rawdata$index$t
	validObject(dt)
	
	return(data.table(mf))
	