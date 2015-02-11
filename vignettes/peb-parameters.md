---
title: "Vignette Title"
author: "Vignette Author"
date: "2014-12-23"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

## Parameters for the Parametric Empirical Bayes Algorighm 

When using the Paremetric Empritical Bayes (PEB) algorithm to set personalized
thresholds for a biomarker, it is necessary estimate the population average
(`mu`) for the biomarker as well as the within (sigma squared, `s2`) and between
(sigma squared, `t2`) person variation of that marker.  Observed marker levels
from a population of healthy individuls (the 'prior population' here on), at
least some of whome have marker levels taken from multiple time points ('prior
observations') is required to estimate each of these paremters.

## Parameter estimates with equal k's 

In the case that each member prior of the population contributes equal numbers
of prior observations, the mean of all prior observations is a minimum variance
unbiased (MVU) estimate of the population mean and the mean of the within person
variances is an MVU estimate of sigma squared, and the `var(m_i) - (s2/k)` is
an MVU estimate of the the `s2'.  Note that we're subtracting
an estimate of the excess variation attributable to the process of extimateing
the within person means.

To illustrate the processe of estimating these paraemters, we can generate a
a random prior population as follows.

'''{r}
# number of individuals in the prior population
N = 20 

# numbers of prior observations per individual
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
'''

The function `PEBparams` calculates the parameters `mu`, `s2`, and `t2` in
exactly this way when the `'iterative'` is specified.

'''{r}
parmas <- PEBparams(x,id,method='iterative',iterations=1)

# The individual parmeters can be extracted via:
params$mu
params$s2
params$t2
'''

Note that in addition to estimating each of the parameters, `PEBparams` also
returned estimates of the error variance for each parameter (`params$t2.ev`,
etc.), which it uses to calculate confidence intervals for each parameter.

## Parameter estimates with un-equal k's 

In practice, the pior population will typically be a convience populatoin, which
is unlikely to be structured with equal numbers of observations from individual.
A convenience population can be simulated as follows:

'''{r}

# number of individuals in the prior population
N = 20 

# distribution of numbers of observations per individual
d = 2:5 

# numbers of prior observations from each individual
k <- sample(d,N,TRUE)

# a vector of person id's
id = factor(rep(seq(N),times=k), levels=seq(N))

# random personal means
x1 <- rep(rnorm(N),times=k)

# random within person variation
x2 <- rnorm(length(id)) 

# the interclass correlation coefficient
beta <- .6

# random marker levels with ICC = beta
x <- sqrt(beta)*x1 + sqrt(1-beta)*x2

'''

**Note that** `d` begins with 2 in order to example to keep the example as simple
as possible.  In principal including individuals with just one prior observation
is possible, and including such individuals improves estimates of `mu` and `t2`. 

With unequal numbers of observations per person, an MVU estimate of `s2` can be
still obtained by as a weighted average of the within person variances, where
the weights are determined by the number of obsersvations from that individual,
as follows:

'''{r}

# s2 estimates for the each individual
s2_i 	<- tapply(x,id,var)

# otpimal weights
w	<- (k - 1) / 2

# an MVU estimate of s2
s2 <- sum(w*s2_i) / sum(w)

'''

and this is is exatly the estimate of `s2` that is produced by `PEBparams()' in
via the iterative method.

'''{r}
params <- PEBparams(x,id,method='iterative')
params$s2
'''

Unfortunatly, closed form MVU estimates of the `mu` and `s2` are not available.
The `mu` parameter is most accuratly thought of as the average of the personal
averages. With unequal numbers of prior observations, an MVU estimat of `mu`
will be a weighted average of the individual's personal means, as in:

'''{r}
mu_i <- tapply(x,id,mean)
'''

The relative weights applied to each personal mean would depend on numbers or
prior observations per person as well as relative magnituds of `t2` and `s2`.
Hence an MVU estimate of `mu` must rely on an MVU estimate of `t2`.  On the
other hand, the estimate of `t2` depends on an estimate of `mu` around which the
sum of the squares of `mu_i` is calculated. 

Hence, the iterative prcoess begins with a naive estiate of `mu` and `t2`,
and iteravely rew-calculates one pareter and then the other to approach an MVU
estimate of each assymtotically. Details of each calculation can be found here
(TODO: make 'here' a link to the PDF the details...).

Note that in the case of equal numbers of observations from each individual in
the prior population, subsequent iterations of the estimation process do not
change the parameter estimates, and hence, in the first example, we could have
specified `PEBparams(...,iterations=10)` with no effect.

## Application of the PEB algorithm

When using the PEB to assess whether or not an individuals marker level is
higher (or lower) than expected given her prior marker levels and the
distribution of markers in the prior population, a threshold (quantile) is
calculated such such that the probability that her current marker level exceeds
that threshold is `p = (1 - specificity)`.  The data required to individual's
threshold depends on the estimated parameters (`mu`,`s2`, and `t2`), the number
of prior observations from that individual (`n`) and the average (mean) of her
prior observations (`ybar`). 

To continue the previous example, 

'''{r}
# number of random markers levels to generate for our ficticous person 
y_k <- 6

# randomly generate a 'true' personal mean 
y_mean <- rnorm(1)*sqrt(beta)

# residuals around her 'true' personal mean 
y_resid <- rnorm(y_k)*sqrt(1-beta)

# a vector of random marker levels for an individual from our ficticious
# population
y <- y_mean + y_resid

# divide the vector into prior and current marker levels
(y_current <- y[y_k])
(y_past <- y[-y_k])

'''

Now, in order to evaluate whether or not the `y_current` is unexpectedly high,
with specficity 95%, we can set a set a threshold, using `qpeb()`, and test
whether or not `y_current` exceeds the resulting threshold, which is caluated
via:

'''{r}
threshold <- qpeb(p=0.95,# set threshold with 95% specificity
			n=length(y_past),
			ybar=mean(y_past),
			mu=params$mu,
			sigma2=params$t2,
			tau2=params$s2)

currentIsElevated <- (y_current > threshold)
'''

If the variable `currentIsElevated` is `TRUE`, then  y_current is considered
'elevated' with 95% specificity by the PEB.

Note that we could have simply passed the `params` object to `qpeb()` instead of
specifying mu,sigma2, and tau2 explicitly, as in,


'''{r}
threshold <- qpeb(p=0.95,# set threshold with 95% specificity
			n=length(y_past),
			ybar=mean(y_past),
			params)

'''

or we could have calculated the probability of the current value using `ppeb()`
as in:

'''{r}
prob <- ppeb(y_current,
			n=length(y_past),
			ybar=mean(y_past),
			params)

currentIsElevated <- (prob > 0.95)
'''

Easeir still, we cold have calculated probablities for the entire series of
marker levels for our ficticious person via: 

'''{r}
probs <- ppeb(y,params)
currentIsElevated <- (probs[y_k] > 0.95)
'''

Notice that if the paramters `n` and `ybar` are not provided, `ppeb` calculates
`n` and `ybar` for each value in y using the assumption that y is a vector containing a
series of marker results from a **single** individual.


## Sample Size requirements

In practice, reasons for setting thresholds for an assessment include (1) that
insufficiently elevates marker levels are not indicitive of a given outcome
(disease), and (2) that a positive test may require follow up exaiminations or
treatments which involve risks that an individual who is not at sufficient risk
of the outcome of interest will want to avoid.  For that reason, it is necessary
to ensure that the range of False Positive Rates (FPRs) that are likely to
result from the application of the PEB algorithm is reasonably narrow.  

Devaitions from the expected FPR (=1 - specificity) may result from errors in
the estimates of the parmeters `mu`, `s2`, and `mu`, as well as from failure of
the biomarker to be distributes normally within individuals, between,
individuals, or both.  Failure of the biomarker meet assumtptions on it's
distribution is beyond the topic of this discussion. 

Because the effective sample sizes used to estimate the `t2` and `s2` parameters
are smaller than that of the `mu` parameter, errors in `t2` and `s2` are likeley
to make greater contributions to deviations from the expected FPR.  (See the PDF
for discussion on the effective sample sizes of each parameter estimate.)
Fortunately, the function `qpeb` can provide confidence intervals around the
expected FPR.  **Note that the confidence intervals provided by `qpeb` account
only for deviations in the estiamtes of the parameters `t2` and `s2` and not for
the deviations in the estimate of `mu`.**


'''{r}
qpeb(p = 0.95,
	params,
	n=seq(0,y_k-1),
	ybar=c(0,cumsum(y)/seq(y_k))[-y_k + 1],
	conf.level=0.9)
'''

Note that the confidence intervals on the FPR depend on (1) the size and
structure of the prior population, (2) the ICC, and (3) the length of the
history.  Notably, the confidence intervals are not dependent on the
individual's personal average, so we could have let `ybar=0` in the above call
to `qpeb`

## Comparison of the 'iterative' and 'ICC1' methods

The parameters `s2` and `t2` are related to the Interclass Correlation
Coeffiect (beta) by the equaion 'beta = t2/(s2 + t2)', which can be estimated
using ANOVA methods via





