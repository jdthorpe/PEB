## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----,echo=FALSE---------------------------------------------------------
#-- obs <- pmax(rpois(20,lambda=2.5),1)
#-- within <- rnorm(sum(obs)) 
#-- between <- rep(rnorm(20),obs)
#-- person = c(1,1,2,3,3,3,4,4,4,5)
#-- myData <- data.frame(person = factor(rep(letters[1:20],obs)),
#-- 					markerLevels = round(100 + within *3+ between*7,1))
#-- dput(myData)

myData <- structure(list(person = structure(c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 
3L, 3L, 4L, 4L, 5L, 5L, 5L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 8L, 8L, 
8L, 8L, 8L, 8L, 9L, 9L, 10L, 11L, 11L, 12L, 12L, 13L, 13L, 13L, 
13L, 13L, 13L, 14L, 14L, 15L, 15L, 16L, 16L, 17L, 17L, 18L, 18L, 
18L, 18L, 19L, 19L, 20L, 20L, 20L), .Label = c("a", "b", "c", 
"d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", 
"q", "r", "s", "t"), class = "factor"), markerLevels = c(104.6, 
98.8, 109.7, 107.7, 112.9, 102.7, 106.7, 110.4, 105.4, 93.8, 
91.3, 100.2, 105, 102, 100.5, 101.7, 103.6, 101.9, 102.3, 104.5, 
108.4, 100.6, 102.8, 97, 103.1, 100.1, 101.5, 102.6, 103.6, 104.2, 
99.4, 100.2, 101.2, 99, 104.8, 100.5, 104.9, 103.9, 102.4, 103.9, 
88, 87.8, 99.5, 91.7, 108.5, 103.1, 100.9, 104.3, 100.5, 95.5, 
96.4, 94.9, 99.2, 99.6, 93.9, 93.4, 94.8)), .Names = c("person", 
"markerLevels"), row.names = c(NA, -57L), class = "data.frame")

## ------------------------------------------------------------------------
library('PEB')
head(myData,10)
params <- PEBparams(x=myData$markerLevels,
					id=myData$person)
params

## ------------------------------------------------------------------------
# let's use person 'a' as an example.  Her prior values are:
(y_past = myData$markerLevels [myData$person == 'a'])

# and a personalized threshold with 95% specificity can be set via:
qpeb(p=0.95,# specificity level
	 n=length(y_past),
	 ybar=mean(y_past),
	 mu=params$mu,
	 sigma2=params$s2,
	 tau2=params$t2)

# or simply:
qpeb(p=0.95,# specificity level
	 n=length(y_past),
	 ybar=mean(y_past),
	 params)

## ------------------------------------------------------------------------
y_current <- 105.8
ppeb(y_current,
	 n=length(y_past),
	 ybar=mean(y_past),
	 params)

## ------------------------------------------------------------------------
y <- c(y_past,
	   y_current)
ppeb(y,params)

## ------------------------------------------------------------------------
# calculate the parameters using the Iterative method
params_iterative <- 
	PEBparams(x=myData$markerLevels,
			  id=myData$person,
			  method='iterative')

qpeb(p = 0.95,
	params_iterative,
	n=length(y_past),
	ybar=mean(y_past),
	conf.level=0.9)

