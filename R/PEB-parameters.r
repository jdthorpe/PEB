# --------------------------------------------------------------------------------------
# Programmer: Jason Thorpe
# Date        03/17/2011
# Language:   R (Version 2.12.0)
# Purpose:    
# Comments:   
# --------------------------------------------------------------------------------------

#' Estimate the PEB parameters
#'
#' Estimate the PEB parameters from a set of marker values (x) from individuals
#'
#' @export
#' @param x observations of the for the variable of intereset
#'
#' @param id a 'factor' with identifiers for the individuals who contributed 
#' the observations in \code{x}
#'
#' @param method Either 'ICC1' or 'iterative' (matched via match.arg). 
#' If method == 'ICC1' (default) multilevel::ICC1 is used to estimate the 
#' within and between group variation. Otherwise an iterative approach is used.
#'
#' @param iterations The number of iterations for estimating the tau 
#' squared parameter (and it's error variance) (default = 2)
#'
#' @param conf.level The width of the confidence interval for sigma-squared 
#' and tau-squared  
#'
#' @family peb
#' @return a list with the following values:
#' s2: sigma quared
#' s2.ev: the error variance for sigma squared
#' t2: tau squared
#' t2.ev: the error variance for tau squared
#' mu: mu 
#' mu.ev: the error variance for mu
#' m: the number of individuals contributing to the estimates
#' N: the nubmer so observations contributing to the estimates
#' mu.ci: confidence interval for mu
#' t2.ci: confidence interval for tau squared
#' s2.ci: confidence interval for sigma squared
#' nstar: an average sample size for use in calcualing the CI for tau sauqred
#'
#' @examples
#' n = 10
#' m = 200
#' PEBparams(x = rnorm(m*n) + rep(rnorm(m),each=n), id=gl(m,n,m*n))
#'
#' @importFrom multilevel ICC1

PEBparams <- function(x,
					  id,
					  method='ICC1',
					  iterations =2 ,
					  conf.level = 0.9,
					  verbose=FALSE){

	method = match.arg(method,c('ICC1','iterative'))
	if(method == 'ICC1' && !missing(iterations))
		stop("parameter 'iterations' is not used when method == 'ICC1'")

	if(!inherits(id,'factor'))
		stop("parameter 'id' must be a factor")
	if(!  all(levels(id) %in% id)){
		warning(paste("unused factor levels in 'id':",
					  paste(levels(id)[!levels(id) %in% id],
							collapse = ', ')))
		id <- factor(as.character(id))
	}


	N <- length(id)
	m <- length(levels(id))

	if(method == 'ICC1'){

		AOV = aov(x~id)
		(B1 = multilevel::ICC1(AOV))
		V  <-  var(x)
		out  <-  list(s2 = (1-B1) * V,
						   t2 = B1 * V,
						   mu = mean(x),
						   m = m,
						   N = N,
						   method=method)
		class(out) <- c('PEBparams','list')
		return(out)

	}else{#(method == )
		# --------------------------------------------------
		# part 1: non-iterative calculations
		# --------------------------------------------------

		#calcualte various constants	
		nresults <- tapply(id,id,length) # results per person
		ybari <- tapply(x,id,mean)

		#############################################
		# estiamte sigma^2
		#############################################
		s2i 	<- tapply(x,id,var)
		stopifnot(is.na(s2i) == (nresults == 1))

		#relative error variances of 
		s2i_weights	<- (nresults - 1) / 2
		stopifnot( (s2i_weights ==0) == is.na(s2i))

		s2 <- sum((s2i_weights[!is.na(s2i)])*(s2i[!is.na(s2i)])) / sum(s2i_weights[!is.na(s2i)])
		s2ev <- 2*(s2^2)/(N - m)

		rm(s2i,s2i_weights)

		# --------------------------------------------------
		# --------------------------------------------------
		# part 2: iterative calculations of mu and Tau^2
		# --------------------------------------------------
		# --------------------------------------------------

		# ------------------------------
		# part 2a: fixed parameters for the iterative processes
		# ------------------------------
		tr <- table(nresults)
		uk <- as.numeric(names(tr[tr>1]))
		rm(tr)

		# estimates of the t2k's 
		t2k   <- numeric(0)
		for(k in uk)
			t2k 	<- c(t2k, var(ybari[nresults == k] ) - (s2/k) )

		# initialize the variables involved in the iterative process
		t2kev_list  <-  list()
		t2ev  <-  muev  <-  t2  <-  numeric(0)

		# ------------------------------
		# part 2b: initial estiamte of EV(T_k^2) from the tk's
		# ------------------------------
		# Next, estimate the error variances which DOES depend on the estimate for tau^2
		# We use the unweighted average of the t2k's because the t2k's are unbaiased.
		# We cannot weight them according to their value because then underestimates 
		# become overweitghted and underestimated become underweighted and the result 
		# is a biased estimate of t2k and of it's error variance. we do however weight 
		# them according to their m_k's becaue there is not bias there.

		t2kev_temp <- numeric(0)#initialization
		for(k in uk){# this would be more efficient with an apply stmt.
			m_k 		<- sum(nresults == k)
			#in the alternate version we estimate all the t2k's then estimate the t2kev's
			#this estimate DOES depend on the estimate for tau^2
			t2kev_temp 	  <- c(t2kev_temp, 
									2*((mean(t2k) + s2/k)^2) / (m_k - 1) + (s2ev/(k^2))
									)
		}
		t2kev_list  <-  list(t2kev_temp)
		rm(k,m_k)

		# ------------------------------
		# part 2c: iterative re-estimation of Tau2 and mu
		# ------------------------------
		t2 <- mu <- t2ev <- muev <- numeric(0)
		# This function nested for lexical scoping 
		t2_mu  <-  function(){
			#re-estimate t2 using the latest t2kev estimates
			t2new <- sum((1/(t2kev_list[[j]]))*(t2k))/sum(1/(t2kev_list[[j]]))

			#re-estimate t2_ev using the latest t2 estimates
			t2kev_temp <- numeric(0) #initialize
			for(k in uk){
				m_k 		<- sum(nresults == k)
				t2kev_temp 	<- c(t2kev_temp, 
									2*((t2new + s2/k)^2) / (m_k - 1) + (s2ev/(k^2)))
			}
			t2kev_list <<- c(t2kev_list,list(t2kev_temp))

			#############################################
			# estimate mu (which depends on the estimate for tau^2)
			#############################################
			ybari_weights	<- 1/(t2new + (s2/nresults))
			bn <- (t2new)/((t2new) + (s2/nresults))

			# extend the lists of estmates
			t2 <<- c(t2,t2new )
			mu <<- c(mu, sum(ybari_weights*ybari) / sum(ybari_weights))

			t2ev <<- c(t2ev, 1/sum(1/t2kev_temp)) #estimate t2ev from the new t2kev's
			muev <<- c(muev, (t2new)/sum(bn))

		}

		stopifnot(iterations>0)

		probs <- c(((1-conf.level)/2),1-((1-conf.level)/2)) # probabilities for the CI

		for(j in 1:iterations){
			t2_mu()
			if(verbose){
				nstar   <- ((m-length(uk))/sqrt((t2ev[j]) / (2*((t2[j])^2)/(m-length(uk)))))
				cat('iteration:',j,
					'| Tau Squared CI width: ',
					format(diff((t2[j]*(nstar))/qchisq(p=probs[2:1],df = nstar)),8),
					'|  Mu CI width: ',
					format(diff(mu[j] + sqrt(muev[j])*qnorm(p=probs)),8),'\n')
			}
		}

		#############################################
		# return the estimates
		#############################################

		nstar   <- ((m-length(uk))/sqrt((t2ev[j]) / (2*((t2[j])^2)/(m-length(uk)))))

		out <- list(s2 = s2,
				   s2.ev = s2ev,
				   t2 = t2[iterations],
				   t2.ev = t2ev[iterations],
				   mu = mu[iterations],
				   mu.ev = muev[iterations],
				   m = m,
				   N = N,
				   mu.ci = mu[iterations] + sqrt(muev[iterations])*qnorm(p=probs),
				   t2.ci = t2[iterations]*qchisq(p=probs[2:1],df = nstar)/nstar,
				   s2.ci = s2*qchisq(p=probs[2:1],df = N-m  )/(N-m),
				   nstar = nstar,
				   conf.level = conf.level, 
				   method=method)
		class(out) <- c('PEBparams','list')
		return(out)

	}

}


#' @export
print.PEBparams <- function(x,digits=4){
	cat('mu:',format(x$mu[1],digits=digits))
	if(x$method == 'iterative')
		cat(' [',format(x$mu.ci[1],digits=digits),
		   	',',format(x$mu.ci[2],digits=digits),']')
	cat('\n')

	cat('sigma^2:',format(x$s2[1],digits=digits))
	if(x$method == 'iterative')
		cat(' [',format(x$s2.ci[1],digits=digits),
		   	',',format(x$s2.ci[2],digits=digits),']')
	cat('\n')

	cat('tau^2:',format(x$t2[1],digits=digits))
	if(x$method == 'iterative')
		cat(' [',format(x$t2.ci[1],digits=digits),
		   	',',format(x$t2.ci[2],digits=digits),']')
	cat('\n')

	cat('Calculated using',x$N,'observations, from ',x$m, 'inidividuals.\n',sep = ' ')
	if(x$method == 'iterative')
		cat('([x,y] indicates the ',format(100*x$conf.level,2),'% confidence interval for each parameter)\n',sep = '')
}


qpeb.alt<-function(p,...,n,ybar,conf.level,details=FALSE)  {
	with(pebVarArgs(...),{
				bn <- bn(n)
				sd_n <<- sqrt(sigma + (tau*(1-bn)))
				mu_n <<- mu*(1-bn) + ybar*bn
				})
	if(length(p) == 1)
		return(qnorm(p,mean = mu_n,sd=sd_n))
	out <- matrix(NA,length(n),length(p))
	for(i in 1:length(p))
		out[,i] <- qnorm(p[i],mean = mu_n,sd=sd_n)
	dimnames(out) <- list(if(!is.null(names(n)))names(N) else seq(length(n)),
						  paste('p =',p))
	return(out)
}

