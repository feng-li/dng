##' @title the quantile function of split-normal
##' @param p vector of probabilities.
##' @param mu vector of means.
##' @param sigma vector of standard deviations.
##' @param lmd vector of skewness.
##' @return vector of quantile

qsplitn <- function(p,mu,sigma,lmd)
{
  #quantile = rep(NA,length(p))
  if(p<=(1/(1+lmd)))
  {
    p0 = (1+lmd)*p/2
    quantile = qnorm(p0,mu,sigma)
  }
  if(p>(1/(1+lmd)))
  {
    p0 = (p-(1-lmd)/(1+lmd))*(1+lmd)/(2*lmd)
    quantile = qnorm(p0,mu,(sigma*lmd))
  }

  return(quantile)

}


