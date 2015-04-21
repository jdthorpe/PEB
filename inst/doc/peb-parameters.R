## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
# The number of individuals in the prior population
N = 20 

# The numbers of prior observations per individual
k = 5 

# a vector of person id's
id = factor(rep(seq(N),times=k), levels=seq(N))

# random personal means
x1 <- rep(rnorm(N),times=k)

# random within person variation
x2 <- rnorm(N*k) 

# the interclass correlation coefficient, defined as 't2/(s2 + t2)'
beta <- .6

# random marker levels with ICC = beta
x <- sqrt(beta)*x1 + sqrt(1-beta)*x2

mu = mean(x)
s2 = mean(tapply(x,id,var))
t2 = var(tapply(x,id,mean))- (s2/k)
c(mu,s2,t2)

## ------------------------------------------------------------------------
library('PEB')
params <- PEBparams(x,id,method='iterative',iterations=1)
params

# The individual parameters can be extracted via:
c(params$mu,params$s2,params$t2)

## ------------------------------------------------------------------------
params[c('mu.ev','t2.ev','s2.ev')]

## ------------------------------------------------------------------------
# number of individuals in the prior population
N = 20 

# distribution of numbers of observations per individual
d = 2:5 

# numbers of prior observations from each individual
k <- sample(d,N,replace=TRUE)

# a vector of person id's
id = factor(rep(seq(N),times=k), levels=seq(N))

# random personal means
x1 <- rep(rnorm(N),times=k)

# random within person variation
x2 <- rnorm(sum(k)) 

# the interclass correlation coefficient
beta <- .6

# random marker levels with ICC = beta
x <- sqrt(beta)*x1 + sqrt(1-beta)*x2

## ------------------------------------------------------------------------
# s2 estimates for the each individual
s2_i 	<- tapply(x,id,var)

# optimal weights
w	<- (k - 1) / 2

# an MVU estimate of s2
(s2 <- sum(w*s2_i) / sum(w))

## ------------------------------------------------------------------------
params <- PEBparams(x,id,method='iterative')
params$s2

## ------------------------------------------------------------------------
mu_i <- tapply(x,id,mean)

## ------------------------------------------------------------------------
# number of random markers levels to generate for our fictitious person 
y_k <- 6

# randomly generate a 'true' personal mean 
y_mean <- rnorm(1)*sqrt(beta)

# residuals around her 'true' personal mean 
y_resid <- rnorm(y_k)*sqrt(1-beta)

# a vector of random marker levels for an individual from our fictitious
# population
y <- y_mean + y_resid

# divide the vector into prior and current marker levels
(y_current <- y[y_k])
(y_past <- y[-y_k])

## ------------------------------------------------------------------------
threshold <- qpeb(p=0.95,# set threshold with 95% specificity
			n=length(y_past),
			ybar=mean(y_past),
			mu=params$mu,
			sigma2=params$s2,
			tau2=params$t2)

qpeb(p=0.95, n=5, ybar=.2, mu=0, sigma2=.5, tau2=.5)

(currentIsElevated <- (y_current > threshold))

## ------------------------------------------------------------------------
threshold <- qpeb(p=0.95,# set threshold with 95% specificity
			n=length(y_past),
			ybar=mean(y_past),
			params)

## ------------------------------------------------------------------------
prob <- ppeb(y_current,
			n=length(y_past),
			ybar=mean(y_past),
			params)

specificity <- 0.95
(currentIsElevated <- (prob > specificity))

## ------------------------------------------------------------------------
probs <- ppeb(y,params)
(currentIsElevated <- (probs[y_k] > 0.95))

## ------------------------------------------------------------------------
qpeb(p = 0.95,
	params,
	n=seq(0,y_k-1),
	ybar=c(0,cumsum(y)/seq(y_k))[-y_k + 1],
	conf.level=0.9)

## ------------------------------------------------------------------------
PEB.conf.int(p=0.95, n=5, conf.level=0.9, params=params)

## ------------------------------------------------------------------------
library(multilevel)
beta <- multilevel::ICC1(aov(x~id))

## ------------------------------------------------------------------------
V  <-  var(x)
s2 <- (1-beta) * V
t2 <- beta * V
mu <- mean(x)

## ----, eval = FALSE------------------------------------------------------
#  m = 40 # number of individuals in the prior population
#  d = seq(2,4) # distribution of maximum numbers of prior observations per individual
#  n1 = 500 # number of random population structures to simulate
#  n2 = 20 # number of iterations per population structure
#  
#  beta <- c(.3,.5,.7,.9) # ICC parameters
#  
#  # INITIALIZE CONTAINERS
#  out_ICC1 <-
#  out_Iter1 <-
#  out_Iter2 <-
#  out_Iter10 <- NULL
#  
#  flds <- c('s2','t2','mu')
#  
#  for(i in seq(n1)){
#  	# generate the population structure
#  	k <- sample(d,m,replace=TRUE)
#  	.id = rep(seq(m),times=k)
#  	id = factor(.id, levels=seq(m))
#  
#  	for(j in seq(n2)){
#  
#  		tmp_ICC1 <-
#  		tmp_Iter1 <-
#  		tmp_Iter2 <-
#  		tmp_Iter10 <- c()
#  
#  		for(b1 in beta){
#  			#generate the randoms
#  			x1 <- rep(rnorm(m),times=k)
#  			x2 <- rnorm(sum(k))
#  			x <- sqrt(b1)*x1 + sqrt(1-b1)*x2
#  
#  			tmp_ICC1 	<- c(tmp_ICC1   ,unlist(PEBparams(x,id,method='ICC1')[flds]))
#  			tmp_Iter1 	<- c(tmp_Iter1  ,unlist(PEBparams(x,id,method='iterative',iterations=1)[flds]))
#  			tmp_Iter2 	<- c(tmp_Iter2  ,unlist(PEBparams(x,id,method='iterative',iterations=2)[flds]))
#  			tmp_Iter10 	<- c(tmp_Iter10 ,unlist(PEBparams(x,id,method='iterative',iterations=10)[flds]))
#  
#  		}
#  
#  	out_ICC1   <- rbind(out_ICC1   ,tmp_ICC1  )
#  	out_Iter1  <- rbind(out_Iter1  ,tmp_Iter1 )
#  	out_Iter2  <- rbind(out_Iter2  ,tmp_Iter2 )
#  	out_Iter10 <- rbind(out_Iter10 ,tmp_Iter10)
#  
#  	}
#  
#  	if(! (i %% 10))
#  		cat('completed',i*n2,'iterations\n')
#  }
#  
#  # Convenience functions for calculating the SD and bias
#  SD <- function(x,exp)
#  	sqrt(sum((x - exp)^2)/length(x))
#  bias <- function(x,exp){
#  	if(missing(exp))
#  		mean(x)
#  	else
#  		mean(x-exp)
#  }
#  
#  # Summarize the SD's of the s2 estimates by beta
#  V_sd <- s2_sd <- t2_sd  <- mu_sd  <-
#  V_bias <- s2_bias <- t2_bias <- mu_bias   <-  NULL
#  for(i in 1:length(beta)){
#  	s2_bias <- cbind(s2_bias,
#  				  c(ICC1  = bias(out_ICC1  [,3*(i-1)+1],exp = 1-beta[i]),
#  					Iter1 = bias(out_Iter1 [,3*(i-1)+1],exp = 1-beta[i])))
#  	V_bias <- cbind(V_bias,
#  				  c(ICC1  = bias(out_ICC1  [,3*(i-1)+1] + out_ICC1  [,3*(i-1) + 2],exp = 1),
#  					Iter1 = bias(out_Iter1 [,3*(i-1)+1] + out_Iter1 [,3*(i-1) + 2],exp = 1),
#  					Iter2 = bias(out_Iter2 [,3*(i-1)+1] + out_Iter2 [,3*(i-1) + 2],exp = 1),
#  					Iter10= bias(out_Iter10[,3*(i-1)+1] + out_Iter10[,3*(i-1) + 2],exp = 1)))
#  	mu_bias <- cbind(mu_bias,
#  				  c(ICC1  = bias(out_ICC1  [,3*(i-1)+3]),
#  					Iter1 = bias(out_Iter1 [,3*(i-1)+3]),
#  					Iter2 = bias(out_Iter2 [,3*(i-1)+3]),
#  					Iter10= bias(out_Iter10[,3*(i-1)+3])))
#  	s2_sd <- cbind(s2_sd,
#  				  c(ICC1  = SD(out_ICC1  [,3*(i-1)+1],exp = 1-beta[i]),
#  					Iter1 = SD(out_Iter1 [,3*(i-1)+1],exp = 1-beta[i])))
#  	V_sd <- cbind(V_sd,
#  				  c(ICC1  = SD(out_ICC1  [,3*(i-1)+1] + out_ICC1  [,3*(i-1) + 2],exp = 1),
#  					Iter1 = SD(out_Iter1 [,3*(i-1)+1] + out_Iter1 [,3*(i-1) + 2],exp = 1),
#  					Iter2 = SD(out_Iter2 [,3*(i-1)+1] + out_Iter2 [,3*(i-1) + 2],exp = 1),
#  					Iter10= SD(out_Iter10[,3*(i-1)+1] + out_Iter10[,3*(i-1) + 2],exp = 1)))
#  	mu_sd <- cbind(mu_sd,
#  				  c(ICC1  = SD(out_ICC1  [,3*(i-1)+3],exp=0),
#  					Iter1 = SD(out_Iter1 [,3*(i-1)+3],exp=0),
#  					Iter2 = SD(out_Iter2 [,3*(i-1)+3],exp=0),
#  					Iter10= SD(out_Iter10[,3*(i-1)+3],exp=0)))
#  }
#  
#  # set dim names of the output tables.
#  DNs2  <-  list(c('ICC1', 'iterative'),
#  					 paste('beta =',beta))
#  DN  <-  list(c('ICC1', '1 iteration', '2 iterations', '10 iterations'),
#  					 paste('beta =',beta))
#  
#  table1 <- matrix(paste(format(s2_bias,digits=2,scientific=FALSE),
#  								  ' (',format(s2_sd,digits=2),')',
#  								  sep=""),2,4,dimnames=DNs2)
#  table2 <- matrix(paste(format(V_bias ,digits=2),
#  								  ' (',format(V_sd ,digits=2),')'
#  								  ,sep=""),4,4,dimnames=DN)
#  table3 <- matrix(paste(format(mu_bias,digits=2),
#  								  ' (',format(mu_sd,digits=2),')',
#  								  sep=""),4,4,dimnames=DN)

