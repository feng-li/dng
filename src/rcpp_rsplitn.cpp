#include <Rcpp.h>
using namespace Rcpp;
//' random generation for Split-Normal distribution
//'
//' Density, distribution function, quantile function and random generation for
//' the split-normal distribution.
//'
//' The random variable \eqn{y} follows a split-normal distribution, which has
//' density \deqn{F(y) = 1/(1+\lambda)\sigma \sqrt(2/\pi)
//' exp(-(y-\mu)^2)/2\sigma^2 ,if y<=\mu} , \deqn{F(y) = 1/(1+\lambda)\sigma
//' \sqrt(2/\pi) exp(-(y-\mu)^2)/2\sigma ^2 \lambda ^2 ,if y>\mu}
//'
//' where \eqn{\sigma >0} and \eqn{\lambda >0}. The Split-normal distribution
//' reduces to normal distribution when \eqn{\lambda = 1}
//'
//' @param n number of observations. If length(n) > 1, the length is taken to
//' be the number required.
//' @param mu vector of location parameter. (The mode of the density)
//' @param sigma vector of standard deviations.
//' @param lmd vector of skewness parameters (>0). If is 1, reduced to
//' symmetric student t distribution.
//' @return \code{rsplitt} generates random deviates. Invalid arguments will
//' result in return value NaN, with a warning.
//'
//' The length of the result is determined by \code{n} for \code{rsplitn}.
//'
//' The numerical arguments other than n are recycled to the length of the
//' result. Only the first elements of the logical arguments are used.
//' @author Feng Li, Jiayue Zeng
//' @seealso For split-normal distribution, also see: \code{\link{dsplitn}()}
//' \code{\link{psplitn}()} \code{\link{qsplitn}()} \code{\link{split-n}()}
//' @references Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of
//' conditional distributions using smooth mixtures of asymmetric student t
//' densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
//' @keywords distribution asymmetric normal
//' @examples
//'
//'
//' n = 10
//' mu <- c(0,1,2)
//' sigma <- c(0.5,1,2)
//' lmd <- c(1,2,3)
//' rsplitn0 = rsplitn(n,mu,sigma,lmd)
//'
//'
//' @export
// [[Rcpp::export]]
NumericVector rsplitn(int n, NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int i,j;
  NumericVector u(n),out(n);
  for(i = 0; i<n; i++)
  {
    u[i] = R::runif(0,1);
  }

  int a[4];
  a[0] = u.size();
  a[1] = mu.size();
  a[2] = sigma.size();
  a[3] = lmd.size();

  for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
  for(j = a[2];j<n;j++) { sigma[j] = sigma[j-a[2]];}
  for(j = a[3];j<n;j++) { lmd[j] = lmd[j-a[3]];}

  NumericVector mu_long(n),df_long(n),phi_long(n),lmd_long(n);
  NumericVector I0(n),I(n);
  NumericVector p0std(n), y0std(n);



  for(i = 0;i<n;i++)
  {
    I0[i] = (u[i]<=(1/(1+lmd[i])));

    if(I0[i])
    {
      p0std[i] = u[i]*(1+lmd[i])/2;
      y0std[i] = R::qnorm5(p0std[i],mu[i], sigma[i], TRUE, FALSE);
      out[i] = y0std[i]*sigma[i]+mu[i] ;
    }

    else
    {
      p0std[i] = (u[i]-1/(1+lmd[i]))*(1+lmd[i])/(2*lmd[i])+0.5;
      y0std[i] = R::qnorm5(p0std[i], mu[i],sigma[i], TRUE, FALSE);
      out[i] = y0std[i]*(sigma[i]*lmd[i])+mu[i] ;
    }

  }


  return out;
}
