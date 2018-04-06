#include <Rcpp.h>
using namespace Rcpp;
//' quantile function of split-normal distribution
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
//' @param p vector of probability.
//' @param mu vector of location parameter. (The mode of the density)
//' @param sigma vector of standard deviations.
//' @param lmd vector of skewness parameters (>0). If is 1, reduced to
//' symmetric student t distribution.
//' @return Gives the quantile functio Asymmetric Normal distribution. The
//' random variable y follows a split-normal distribution, y~N(\eqn{\mu},
//' \eqn{\sigma}, \eqn{\lmd}).
//' @author Feng Li, Jiayue Zeng
//' @seealso For split-normal distribution, also see: \code{\link{dsplitn}()}
//' \code{\link{psplitn}()} \code{\link{rsplitn}()} \code{\link{split-n}()}
//' @references Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of
//' conditional distributions using smooth mixtures of asymmetric student t
//' densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
//' @keywords distribution asymmetric normal
//' @examples
//'
//'
//' p <- c(0.25,0.5,0.75)
//' mu <- c(0,1,2)
//' sigma <- c(0.5,1,2)
//' lmd <- c(1,2,3)
//' qsplitn0 = qsplitn(p,mu,sigma,lmd)
//'
//'
//' @export
// [[Rcpp::export]]
NumericVector qsplitn(NumericVector p,NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int a[4];
  int n,i,j;
  a[0] = p.size();
  a[1] = mu.size();
  a[2] = sigma.size();
  a[3] = lmd.size();


  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] ) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=3;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { p[j] = p[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { sigma[j] = sigma[j-a[2]];}
    for(j = a[3];j<n;j++) { lmd[j] = lmd[j-a[3]];}
  }

  NumericVector p0(n),quantile(n);

  for(int i=0;i<n;i++)
  {
    if(p[i]<=(1/(1+lmd[i])))
    {
      p0[i] = (1+lmd[i])*p[i]/2;
      quantile[i] = R::pnorm5(p0[i],mu[i],sigma[i],1,0);

    }

    else
    {
      p0[i] = (p[i]-(1-lmd[i])/(1+lmd[i]))*(1+lmd[i])/(2*lmd[i]);
      quantile[i] = R::pnorm5(p0[i],mu[i],(sigma[i]*lmd[i]),1,0);
    }

  }

  return quantile;

}
