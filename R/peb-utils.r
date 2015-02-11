#' Calculate each of the variability estimates used by the PEB
#' 
#' A utility function used to calculate each of the variability 
#' estimates used by the PEB algorithm from from any two of them 
#' 
#' @export
#' @family peb
#' 
#' @param mu The population mean
#' 
#' @param sigma2 The within group variation
#' 
#' @param tau2 The between group variation
#' 
#' @param v The total variation (=sigma2 + tau2)
#' 
#' @param beta1,beta The proportion of toal variation attributable to between group variation (=tau2 / v)
#' 
#' @note mu,sigma2,tau2,v, and beta may be substitured by a the object returned by \code{\link{PEBparams}}
#' 
#' @return a list containing the values sigma2, tau2, v, beta1, and mu, 
#' and a function 'bn' which calculates bn for various values of n 
#' given these parameters in this 
#' 

pebVarArgs <-  function(...){
	# calculates all the 
	ARGS = list(...)
	# fastPath out for the special case that the only argument
	# was the product of PEBparams()
	if((length(ARGS) == 1) & inherits(ARGS[[1]],'PEBparams')){
		params <- ARGS[[1]]
		sigma2 = params$s2
		tau2 = params$t2
		v = sigma2 + tau2
		mu = params$mu
		fun = function(n)tau2/((sigma2/n) + tau2)
		return(list(sigma2=sigma2,
					tau2=tau2,
					v=v,
					mu=mu,
					beta1=fun(1),
					bn = fun))
	}

	sigma2 = if('sigma2' %in% names(ARGS)) ARGS$sigma2 else NA
	tau2 = if('tau2' %in% names(ARGS)) ARGS$tau2 else NA
	v = if('v' %in% names(ARGS)) ARGS$v else NA
	beta1 = (if('beta1' %in% names(ARGS))
				ARGS$beta1 
			else if('beta' %in% names(ARGS))
				ARGS$beta 
			else NA)

	if(!is.na(sigma2)) 	stopifnot(sigma2>0)
	if(!is.na(tau2)) 	stopifnot(tau2>0)
	if(!is.na(v)) 		stopifnot(v>0)
	if(!is.na(beta1)) 	stopifnot(beta1>0)

	if(sum(c((is.na(sigma2)||is.nan(sigma2)),
			 (is.na((tau2)||is.nan(tau2))),
			 (is.na((v)||is.nan(v))),
			 (is.na((beta1)||is.nan(beta1))))) != 2)
		stop("please specifiy exacty two parameters out of 'sigma2','tau2','v', and 'beta'")
	# beta1 / v specified
	if(is.na(tau2) & (!is.na(beta1))& (!is.na(v)))
		tau2 = beta1*v
	if(is.na(sigma2) & (!is.na(beta1))& (!is.na(v)))
		sigma2 = (1-beta1)*v
	# sigma2 /v specified
	if(is.na(sigma2) & (!is.na(tau2))& (!is.na(v)))
		sigma2 = v - tau2
	# tau2/v specified
	if(is.na(tau2) & (!is.na(sigma2))& (!is.na(v)))
		tau2 = v - sigma2
	# tau2 / V specified
	if(is.na(sigma2) & (!is.na(tau2))& (!is.na(beta1)))
		sigma2 = (tau2/beta1) - tau2
	# tau2 / V specified
	if(is.na(tau2) & (!is.na(sigma2))& (!is.na(beta1)))
		tau2 = (sigma2*beta1)/(1-beta1)
	if(is.na(sigma2) |is.na(tau2) )
		stop('Insufficient arguments:  ')
	if(is.na(v))
		v = sigma2 + tau2
	if(is.na(beta1))
		beta1 = tau2 / v
	return(list(sigma2=sigma2,
			   	tau2=tau2,
			   	v=v,
			   	beta1=beta1,
				mu=ARGS$mu,
				bn = function(n)tau2/((sigma2/n) + tau2)))
}

