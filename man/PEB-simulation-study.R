# a simulation sudy of the differeces between the 
# estimates of the mean, total variance and within 
# person variance using either the ICC1 function 
# or the Thorpe method.

# with large sample sizes, both methods should approach 
# perfect accuracy and precision.  Hence I'm using 
# relatively small populations (20 individuals) 
# and relatively few sample per person (between 2 and 6)

m = 20 # number of individuals in the population
n = seq(2,6) # distribution of maximum numbers of observations per individual
n1 = 500 # number of random population structures
n2 = 20 # number of iterations per population structure


beta <- c(.3,.5,.7) # b1 parameters

# INITIALIZE CONTAINERS
out_ICC1 <- 
out_Iter2 <- 
out_Iter5 <- 
out_Iter10 <- NULL

flds <- c('s2','t2','mu')

for(i in seq(n1)){
	# generate the population structure
	samplesPerPerson <- sample(n,m,TRUE)
	.id = rep(seq(m),times=samplesPerPerson)
	id = factor(.id, levels=seq(m))

	for(j in seq(n2)){

		tmp_ICC1 <- 
		tmp_Iter2 <- 
		tmp_Iter5 <- 
		tmp_Iter10 <- c()

		for(b1 in beta){ 
			#generate the randoms
			x1 <- rep(rnorm(m),times=samplesPerPerson)
			x2 <- rnorm(length(id)) 
			x <- sqrt(b1)*x1 + sqrt(1-b1)*x2

			tmp_ICC1 	<- c(tmp_ICC1   ,unlist(PEBparams(x,id,method='ICC1')[flds]))
			tmp_Iter2 	<- c(tmp_Iter2  ,unlist(PEBparams(x,id,method='It',iterations=2)[flds]))
			tmp_Iter5 	<- c(tmp_Iter5  ,unlist(PEBparams(x,id,method='It',iterations=5)[flds]))
			tmp_Iter10 	<- c(tmp_Iter10 ,unlist(PEBparams(x,id,method='It',iterations=10)[flds]))

		}

	out_ICC1   <- rbind(out_ICC1   ,tmp_ICC1  )
	out_Iter2  <- rbind(out_Iter2  ,tmp_Iter2 )
	out_Iter5  <- rbind(out_Iter5  ,tmp_Iter5 )
	out_Iter10 <- rbind(out_Iter10 ,tmp_Iter10)

	}

	if(! (i %% 10))
		cat('completed',i*n2,'iterations\n')
}

apply(out_ICC1  ,2,mean)
apply(out_Iter2 ,2,mean)
apply(out_Iter5 ,2,mean)
apply(out_Iter10,2,mean)

indx  <-  seq(0,length(beta)-1) * length(beta)

s2_percent_bias <- 
	rbind(ICC1  = (apply(out_ICC1  [,indx + 1],2,mean)-(1-beta))/(1-beta),
		  Iter2 = (apply(out_Iter2 [,indx + 1],2,mean)-(1-beta))/(1-beta),
		  Iter5 = (apply(out_Iter5 [,indx + 1],2,mean)-(1-beta))/(1-beta),
		  Iter10= (apply(out_Iter10[,indx + 1],2,mean)-(1-beta))/(1-beta))

V_percent_bias <- 
	rbind(ICC1  = apply(out_ICC1  [,indx + 1] + out_ICC1  [,indx + 2],2,mean)-(1),
		  Iter2 = apply(out_Iter2 [,indx + 1] + out_Iter2 [,indx + 2],2,mean)-(1),
		  Iter5 = apply(out_Iter5 [,indx + 1] + out_Iter5 [,indx + 2],2,mean)-(1),
		  Iter10= apply(out_Iter10[,indx + 1] + out_Iter10[,indx + 2],2,mean)-(1))

mu_bias <- 
	rbind(ICC1  = apply(out_ICC1  [,indx + 3],2,mean),
		  Iter2 = apply(out_Iter2 [,indx + 3],2,mean),
		  Iter5 = apply(out_Iter5 [,indx + 3],2,mean),
		  Iter10= apply(out_Iter10[,indx + 3],2,mean))


# calculat the empirical standarad errors
SE <- function(x,exp)
	sqrt(sum((x - exp)^2)/length(x))

# summarize the SE's of the s2 estimates by beta
V_se <- s2_se <- t2_se <- mu_se <- NULL
for(i in 1:length(beta)){
	s2_se <- cbind(s2_se,
				  c(ICC1  = SE(out_ICC1  [,i - 1 + 1],exp = 1-beta[i]),
					Iter2 = SE(out_Iter2 [,i - 1 + 1],exp = 1-beta[i]),
					Iter5 = SE(out_Iter5 [,i - 1 + 1],exp = 1-beta[i]),
					Iter10= SE(out_Iter10[,i - 1 + 1],exp = 1-beta[i])))

	V_se <- cbind(V_se,
				  c(ICC1  = SE(out_ICC1  [,i - 1 + 1] + out_ICC1  [,i - 1 + 2],exp = 1),
					Iter2 = SE(out_Iter2 [,i - 1 + 1] + out_Iter2 [,i - 1 + 2],exp = 1),
					Iter5 = SE(out_Iter5 [,i - 1 + 1] + out_Iter5 [,i - 1 + 2],exp = 1),
					Iter10= SE(out_Iter10[,i - 1 + 1] + out_Iter10[,i - 1 + 2],exp = 1)))

	mu_se <- cbind(mu_se,
				  c(ICC1  = SE(out_ICC1  [,i - 1 + 3],exp = 0),
					Iter2 = SE(out_Iter2 [,i - 1 + 3],exp = 0),
					Iter5 = SE(out_Iter5 [,i - 1 + 3],exp = 0),
					Iter10= SE(out_Iter10[,i - 1 + 3],exp = 0)))
}

names(s2_percent_bias) <- 
names(V_percent_bias) <- 
names(mu_bias) <- 
names(V_se) <- 
names(mu_se) <- 
names(s2_se) <- list(c('ICC1','2 iterations','5 iterations','10 iterations'),
					 paste('beta =',beta))

# the absolute biases and standard errors
s2_percent_bias
V_percent_bias
mu_bias
mu_se
s2_se
V_se

# CALCULATES THE ABSOLUTE RELATIVE ERROR
relative <- function(x) 
	abs(t(t(x)/x[1,]	))

# THE THORPE METHOD HAS SUBSTANTIALLY LOWER BIAS
relative(s2_percent_bias)
relative(V_percent_bias)
relative(mu_bias)

# BOTH METHODS HAVE SIMILAR STANDARD ERRORS
relative(mu_se)
relative(s2_se)
relative(V_se)




