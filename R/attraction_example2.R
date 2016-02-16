#source('attraction.R')
#source('itersur.R')
#source('attraction_data.R')


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

#source('attraction.R')
##source('itersur.R')
#source('attraction_data.R')
