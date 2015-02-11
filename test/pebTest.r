require(PEB)

local({

	# ------------------------------------------------------------
	# well tested functions taken from
	# r:\Urban_N\Rivkin\NovelMarkers\Common_Files\Programs\Quality_Control\NMT\ PEB\ QC\ functions.r
	# ------------------------------------------------------------

	ppeb2 <- function(x,mu,sig,tau,n,ybar){
	# this peb function calculates the PEB in line, which is to say
	# that it assumes that n and ybar have been calculated and 
	# are appropriate for the current value x
		bn <- tau/((sig/n) + tau)
		yn_hat <- mu*(1-bn) + ybar*bn
		vn_hat <- sig + (tau*(1-bn))
		z <- (x - yn_hat)/sqrt(vn_hat)
		p <- pnorm(z)
		return(p)
	}

	qpeb2 <- function(p,mu,sig,tau,n,ybar){
		# this peb function calculates quantiles for the PEB in line, which is to say
		# that it assumes that n and ybar have been calculated and 
		# are appropriate for the current value quantile
		bn <- tau/((sig/n) + tau)
		yn_hat <- mu*(1-bn) + ybar*bn
		vn_hat <- sig + (tau*(1-bn))	
		return(yn_hat + qnorm(p)*sqrt(vn_hat))
	}

	# ------------------------------------------------------------
	# load the PEB package
	# ------------------------------------------------------------


	nvalidations <- 100
	for(i in 1:nvalidations){
		if(!(i %% ceiling(nvalidations/10)))
		cat('i =',i,'\n')

		x <- rnorm(10)
		(tmp <- PEB::ppeb(x,mu=0,sigma=1,tau=0.5,details=T))

		# validate the PEB p-value calculations
		stopifnot(all(abs(tmp$p - ppeb2(x,
										mu=0,
										sig=1,
										tau=0.5,
										n=tmp$n,
										ybar=tmp$ybar) )< 1E-15)) 

		# validate the PEB threshold calculations
		p <- runif(1)
		PEB::qpeb(p,mu=0,sigma=1,tau=0.5,n=tmp$n,ybar=tmp$ybar)
		stopifnot(all(abs(PEB::qpeb(p,mu=0,sigma=1,tau=0.5,n=tmp$n,ybar=tmp$ybar) - 
						  qpeb2(p,
								mu=0,
								sig=1,
								tau=0.5,
								n=tmp$n,
								ybar=tmp$ybar) )< 1E-15))

	}

})


