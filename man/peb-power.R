# ----------------------------------------
# sample size calculations
# ----------------------------------------
# 'power' for setting a PEB threshold really means the 
# range of true False Positive Rates (FPR) that you 
# can expect when setting a threshold with a given 
# specificity for a woman with a given number of prior
# observations.  For example, if you expecd a very high
# specificity, because the intervention that follows a
# positive result is in vasive, you will want to take 
# steps to ensure that the PEB is unlikely to have a higher
# than expected false positive rate.  
#
# The range of actual FPRs depends on the number of samples
# that are used to calculate parameters for the PEB. Technically
# the parameters for the PEB can be derived from as few as 4 samples
# but the parameters would likely be very poor and could 
# result in a FPR which is *much* higher (or lower) than expected.
#
# The script below illustrates calculation of confidence 
# intervals describing the range of FPRs that can be expected
# when using parameters from observations that are avialiable 
# to you prior to your study.
#
# *NOTE* that the confidence intervals provided by qpeb() 
# represent error attributable to estimation of the within and between 
# group variance parameters (\eqn{sigma^2} and \eqn{tau^2}) the variation that can be expected
# a large simualtion study evaluating the effects of the 
# using the PEB in your specific study is recommended.

m = 20 # number of individuals in the population
n = seq(2) # distribution of numbers of observations per individual
n1 = 5 # number of random population structures
n2 = 20 # number of iterations per population structure

nhist  <-  0:5 # numbers of samples in women't history
conf.level  <-  0.95
prob = 0.95 # 1- FPR

beta <- c(.3,.5,.7) # b1 parameters

.out  <-  NULL

for(i in seq(n1)){
	# generate the population structure
	samplesPerPerson <- sample(n,m,TRUE)
	.id = rep(seq(m),times=samplesPerPerson)
	id = factor(.id, levels=seq(m))

	for(j in seq(n2)){
		tmp  <-  NULL

		for(b1 in beta){ 
			#generate the randoms
			x1 <- rep(rnorm(m),times=samplesPerPerson)
			x2 <- rnorm(length(id)) 
			x <- sqrt(b1)*x1 + sqrt(1-b1)*x2

			params  <-  PEBparams(x,id,method='It',iterations=3)
			
			tmp <- cbind(tmp, 
						 qpeb(p=prob,
							  n = nhist,
							  ybar=0,# this parameter falls out...
							  params, conf.level = conf.level)[,-1])
		}
		paste(rep(paste('b1=',beta),each=2),
			  c('upperCI','lowerCI'))
		.out   <- rbind(.out,
					   cbind(nhist,tmp  ))
	}

	if(! (i %% 10))
		cat('completed',i*n2,'iterations\n')
}

OUT  <-  NULL
for(i in seq(2*length(beta)))
	OUT  <-  cbind(OUT, tapply(out[,1+i], .out[,1], median,na.rm=T))

colnames(OUT)  <-  paste(rep(paste('b1=',beta),each=2), c('upperCI','lowerCI'))
OUT



