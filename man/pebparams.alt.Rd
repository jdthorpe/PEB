\name{PEBparams}
\alias{PEBparams}
\title{given these parameters in this
Estimate the PEB parameters}
\usage{
  PEBparams(x, id, method = "ICC1", iterations = 2,
    conf.level = 0.9, verbose = FALSE)
}
\arguments{
  \item{x}{observations of the for the variable of
  intereset}

  \item{id}{a 'factor' with identifiers for the individuals
  who contributed the observations in \code{x}}

  \item{method}{Either 'ICC1' or 'iterative' (matched via
  match.arg). If method == 'ICC1' (default)
  multilevel::ICC1 is used to estimate the within and
  between group variation. Otherwise an iterative approach
  is used.}

  \item{iterations}{The number of iterations for estimating
  the tau squared parameter (and it's error variance)
  (default = 2)}

  \item{conf.level}{The width of the confidence interval
  for sigma-squared and tau-squared}
}
\value{
  a list with the following values: s2: sigma quared s2.ev:
  the error variance for sigma squared t2: tau squared
  t2.ev: the error variance for tau squared mu: mu mu.ev:
  the error variance for mu m: the number of individuals
  contributing to the estimates N: the nubmer so
  observations contributing to the estimates mu.ci:
  confidence interval for mu t2.ci: confidence interval for
  tau squared s2.ci: confidence interval for sigma squared
  nstar: an average sample size for use in calcualing the
  CI for tau sauqred
}
\description{
  Estimate the PEB parameters from a set of marker values
  (x) from individuals
}
\examples{
n = 10
m = 200
PEBparams(x = rnorm(m*n) + rep(rnorm(m),each=n), id=gl(m,n,m*n))
}
\seealso{
  \code{link{PEB}}
}

