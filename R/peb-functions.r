# --------------------------------------------------------------------------------------
# Programmer: Jason Thorpe
# Date        02/24/2011
# Language:   R (Version 2.6.0)
# Purpose:    This package has functions for application of the PEB algorithm
# Comments:
# --------------------------------------------------------------------------------------

#' Functions that apply the Parametric Empirical Bayes (PEB) algorithm
#'
#' Functions that apply the Parametric Empirical Bayes (PEB) algorithm.
#' 
#' \code{zpeb(x,[n,ybar,]...)} returnes the z-score for a marker value after applying the PEB algorithm
#' 
#' \code{ppeb(x,[n,ybar,]...)} returnes the p-values for a marker value after applying the PEB algorithm
#' 
#' \code{qpeb(p,n,ybar,...,conf.level)} returnes the quantiles for a marker if the PEB is 
#' to algorithm with a specificity \code{p}. In particular, if the normality 
#' assumptions are satisified, the probability that a marker exceeds 
#' \code{threshold = qpeb(p,n,ybar,...)} is \code{p}. 
#' 
#' @name peb
#' @family peb
#' 
#' @param x A vector of marker values to be interpreted by the PEB.  \strong{Note}
#' If n and ybar are not specified, x is interpreted to be a series of observed 
#' marker levels from a single individual
#' 
#' @param n \code{n[i]} is the number of prior observations from the same
#' individual which preceed the i'th observation (\code{x[i]}).  Note that 
#' n and ybar are required paremeters for qpeb and for zpeb and ppeb, must 
#' both be present or both be absent.
#' 
#' @param ybar \code{n[i]} is the mean of prior observations from the same 
#' individual which preceed the i'th observation (\code{x[i]}).  Note that 
#' n and ybar are required paremeters for qpeb and for zpeb and ppeb, must 
#' both be present or both be absent.
#' 
#' @param details If True, a data frame containing resutls from each 
#' sub-calculation is returned,otherwise a vector with the z-scores 
#' or p-values is returned. (Default = FALSE)
#' 
#' @param na.rm If TRUE and the parameters \code{n} and \code{ybar} are not passed, then 
#' mising values in X are ignored when calculationg n and ybar from x. \strong{Note}
#' that if \code{n} and \code{ybar} are not passed and all values following the first missing 
#' value in x will also be missing. 
#' 
#' @param conf.level (optional) confidence level for the False Positve Rate (FPR)
#' for the thresholds provided \code{qpeb()}. If conf.level is specified, \code{length(p)}
#' must be 1 and the PEB parameters must be specified by the object returned by 
#' \code{PEBparams()} using the \code{'Iterative'} method, as in:
#' 
#' \code{params = PEBparams(x,id,method='Iterative',iterations=2);
#' qpeb(p,params,conf.level=0.95)}
#' 
#' @inheritParams pebVarArgs
#' 
#' @note Note that each function requires just 2 of the 4 parameters v,sigma2,tau2, and beta (=beta1)
#' \strong{Note} also that mu,sigma2,tau2,v, and beta may be substitured by a the object returned by \code{\link{PEBparams}}
#' 
#' @examples
#' 
#' n = 10
#' m = 200
#' params <- PEBparams(x = rnorm(m*n) + rep(rnorm(m),each=n), id=gl(m,n,m*n))
#' x <- rnorm(10) + 6
#' plot(seq(length(x)),ppeb(x,params))
#' 
#' # interesting facts about the PEB...
#' b1 = .7
#' 
#' # this should be uniformly distributed by definition
#' hist(ppeb(rnorm(5000)*sqrt(1-b1),mu=0,v=1,beta=b1))
#' 
#' # this should be positive on the outset. 
#' hist(ppeb(rnorm(5)*sqrt(1-b1) + 3,mu=0,v=1,beta=b1)) # 10/sqrt(1-b1) standard deviations above normal!!
#' 
#' # but even this should be uniformly distributed in the long run
#' hist(ppeb(rnorm(5000)*sqrt(1-b1) + 10,mu=0,v=1,beta=b1)) # 10/sqrt(1-b1) standard deviations above normal!!
#' 

NULL


#' @export
#' @rdname peb
#' @usage zpeb(x, ..., [n, ybar,] na.rm = FALSE, details = FALSE)
zpeb<- function(x,...,n,ybar,na.rm=FALSE,details=FALSE) {
	if(missing(n) != missing(ybar))
		stop("parameters 'n' and 'ybar' must be either both present or both absent")

	if(missing(n)){
		if(na.rm){
			x.na <- is.na(x)
			ybar <- n <- rep(NA,length(x))
			if(any(!x.na)){
				# personal means
				ybar[!x.na] <- ybarCalc(x[!x.na])
				# number of prior results
				n[!x.na] <- nCalc(x[!x.na])
			}
		}else{
			# personal means
			ybar <- ybarCalc(x)
			# number of prior results
			n <- nCalc(x)
		}
	}
	
	# DO THE WORK
	out <- with(pebVarArgs(...),{
			bn <- bn(n)
			vn_hat <- sigma2 + (tau2*(1-bn))
			sdn_hat <- sqrt(vn_hat)
			yn_hat <- mu*(1-bn) + ybar*bn
			z <- (x - yn_hat)/sdn_hat
			data.frame(x=x,
				  n=n,
				  ybar=ybar,
				  bn=bn,
				  z=z,
				  mu_n=yn_hat,
				  sd_n=sdn_hat,
				  z=z)
			})
	if(details)
		return(cbind(out 
					,p=pnorm(out[,'z'])
					))
	else
		return(out$z)
}

#' @export
#' @usage ppeb(x, ..., [n, ybar,] na.rm = FALSE, details = FALSE)
#' @rdname peb
ppeb<-function(...,details=FALSE)  {
	# determine the probability of a sereies of marker values
	if(details)
		return(zpeb(...,details = T))
	else
		return(pnorm(zpeb(...,details = F)))
}

#-- "
#-- 
#-- source('r:/Urban_N/UrbanGrp/Data Analysis/JThorpe/RTools/PEB/R/peb-parameters.r')
#-- source('r:/Urban_N/UrbanGrp/Data Analysis/JThorpe/RTools/PEB/R/peb-utils.r')
#-- source('r:/Urban_N/UrbanGrp/Data Analysis/JThorpe/RTools/PEB/R/peb-functions.r')
#-- 
#-- params_iter <-  PEBparams(x,id,method='iterative',iterations=2)
#-- params_icc  <-  PEBparams(x,id,method='ICC')
#-- params
#-- 
#-- pebVarArgs(params_iter)
#-- pebVarArgs(params_icc )
#-- 
#-- qpeb(p=0.95, params, n=3, ybar=.5,conf.level=0.95)
#-- 
#-- "

#' @export
#' @usage qpeb(p, ..., n, ybar, verbose = FALSE)
#' @rdname peb
#' @importFrom CompQuadForm davies
qpeb<-function(p,...,n,ybar,conf.level,conf.method='davies', verbose = FALSE)  {

	pebArgs  <-  pebVarArgs(...)
	with(pebArgs,{
				bn_fun <<- bn(n)
				bn <- bn(n)
				bn <<- bn
				s2 <<- sigma2 
				t2 <<- tau2 
				v_n <<- sigma2 + (tau2*(1-bn))
				sd_n <<- sqrt(sigma2 + (tau2*(1-bn)))
				mu_n <<- mu*(1-bn) + ybar*bn
				})
	if(length(p) == 1){
		threshold  <-  qnorm(p,mean = mu_n,sd=sd_n)
		if(missing(conf.level))
			return(threshold)
		params  <-  list(...)[[1]]
		if((!inherits(params,'PEBparams')) || (params$method != 'iterative'))
			stop("If 'conf.level' specified, qpeb must be supplied with a iteragive PEB params object.  (See notes in help file)")

		UN <-  unique(n)
		conf = PEB.conf.int(p,UN,conf.level,params)
		indx  <- match(n,UN)
 		return(cbind(quantile=threshold,
					  p.lower.CI = conf[,'p.lower.CI'][indx],
					  p.upper.CI = conf[,'p.upper.CI'][indx]))
	}

	if(!missing(conf.level))
		stop("'conf.level' not supported when length(p) > 1")

	out <- matrix(NA,length(n),length(p))
	for(i in 1:length(p))
		out[,i] <- qnorm(p[i],mean = mu_n,sd=sd_n)
	dimnames(out) <- list(if(!is.null(names(n)))names(N) else seq(length(n)),
						  paste('p =',p))
	return(out)
}


# just a couple of utility functions so that the 'if(mising(n))){...}' block in zpeb() is not to ugly 
ybarCalc <- function(x)
	return(if(length(x) > 1)
				c(0,cumsum(x[-length(x)])/seq(length(x)-1))
			else 0
			)

nCalc <- function(x)
	return(if(length(x) > 1)
				c(0,seq(length(x)-1))
			else 0
			)


#' Confidence intervals for personalized PEB thresholds 
#' 
#' Confidence intervals for personalized PEB thresholds.
#' 
#' @export
#' @family peb
#' 
#' @param p probability (specificity) 
#' 
#' @param n Number of previous observations in the individual's history. If \code{length(n) == 1}, n is expanded to \code{seq(0,n)}.
#' 
#' @param conf.level Confidence level of the interval.
#' 
#' @param params A PEBparams object returned by \code{PEBparams(...,method='iterative')}
#' 
#' @examples 
#' PEB.conf.int(p=0.95, n=6, conf.level=0.9, params=params)

PEB.conf.int  <-  function(p,n,conf.level,params){
	if((!inherits(params,'PEBparams')) || (params$method != 'iterative'))
		stop("If 'conf.level' specified, qpeb must be supplied with a iteragive PEB params object.  (See notes in help file)")

	if(length(n) == 1)
		n  <-  seq(0,n)

	with(pebVarArgs(params),
	{
		bn <<- bn <- bn(n)
		sd_n <<- sqrt(sigma2 + (tau2*(1-bn)))
	})

	threshold  <-  qnorm(p,sd=sd_n)
	p.lower.CI <- 
	p.upper.CI <- numeric(0)

	with(params,
	{
		alpha <- 1 - conf.level
		cint <- c(alpha/2,1-alpha/2)

		FUN <- function(x,bn,prob){
			CompQuadForm::davies(x,
				   lambda=c(s2/(N-m),t2*(1-bn)/nstar),
				   h = c(N-m  ,nstar),# degrees of freedom

				   # sigma2 *DOES NOT* ACCOUNT FOR THE STANDARD ERROR IN THE ESTIMATE 
				   # OF THE INDIVIDUALS ESTIMATED MEAN WITH N PRIOR OBSERVATIONS
				   # B/C THE PEB THRESHOLD IS SET AS THE SUM OF THE MEAN + THE SQRT 
				   # OF THE WEIGHTED AVERAGE OF TWO CHI-SQUARED VARIALBLES (S2 AND T2)
				   # WHEREAS, THIS sigma2 PARAMETER IS IS ADDED TO THE WEIGHTED AVERAGE 
				   # *WITHIN* THE SQUARE ROOT...
				   # sigma=(1-bn)*mu.se + bn*(sqrt(s2)/n),

				   lim=10000,
				   acc=0.0001)$Qq-prob
		}

		for(i in 1:length(bn)){
			(INT <- (t2*(1-bn[i])*qchisq(p=cint,df = nstar)/(nstar)) + 
					(s2*qchisq(p=cint,df = N-m  )/(N-m  )))
			sd.upper.ci <<- sqrt(uniroot(FUN,
										 interval = INT,
										 bn=bn[i],
										 prob=1-alpha/2)$root)
			sd.lower.ci <<- sqrt(uniroot(FUN,
										 interval = INT,
										 bn=bn[i],
										 prob=alpha/2)$root)
			p.lower.CI[i] <<- pnorm(threshold[i], sd=sd.lower.ci) 
			p.upper.CI[i] <<- pnorm(threshold[i], sd=sd.upper.ci) 
		}
	})
	return(cbind(n, 
				 p.lower.CI = p.lower.CI, 
				 p.upper.CI = p.upper.CI))
}

