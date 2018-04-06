#include <Rcpp.h>
using namespace Rcpp;
//' split-normal distribution
//'
//' Density function for the split normal distribution.
//'
//' Computing the PDF of the Asymmetric Normal distribution. The random
//' variable y follows a split-normal distribution, y~N(\eqn{\mu},
//' \eqn{\sigma}, \eqn{\lmd}), which has density: \deqn{1/(1+\lambda)\sigma
//' \sqrt(2/\pi) exp{-(y-\mu)*2/2\sigma^2}, if y<=\mu}
//' \deqn{1/(1+\lambda)\sigma \sqrt(2/\pi) exp{-(y-\mu)*2/2\sigma^2 \lambda^2},
//' if y>\mu} where \eqn{\sigma>0} and \eqn{\lambda>0}. The Split-normal
//' distribution reduce to normal distribution when \eqn{\lmd=1}.
//'
//' @param x vector of quantiles.
//' @param mu vector of location parameter. (The mode of the density)
//' @param sigma vector of standard deviations.
//' @param lmd vector of skewness parameters (>0). If is 1, reduced to
//' symmetric student t distribution.
//' @param logarithm logical; if TRUE, probabilities p are given as log(p).
//' @return Computing the PDF of the Asymmetric Normal distribution. The random
//' variable y follows a split-normal distribution, y~N(\eqn{\mu},
//' \eqn{\sigma}, \eqn{\lmd}).
//' @author Feng Li, Jiayue Zeng
//' @seealso \code{\link{psplitn}()} for the split-normal distribution.
//' @references Li, F., Villani, M., & Kohn, R. (2010). Flexible modeling of
//' conditional distributions using smooth mixtures of asymmetric student t
//' densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
//' Villani, M., & Larsson, R. (2006) The Multivariate Split Normal
//' Distribution and Asymmetric Principal Components Analysis. Sveriges
//' Riksbank Working Paper Series, No. 175.
//' @keywords distribution asymmetric normal
//' @examples
//'
//' x <- c(0.25,0.5,0.75)
//' mu <- c(0,1,2)
//' lmd <- c(0.5,1,2)
//' sigma <- c(0.5,1,2)
//' dsplitn0 = dsplitn(x, mu, sigma,lmd, TRUE)
//' @export
// [[Rcpp::export]]
NumericVector dsplitn(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector lmd, bool logarithm)
{
  int a[4];
  int n,i,j;
  a[0] = x.size();
  a[1] = mu.size();
  a[2] = sigma.size();
  a[3] = lmd.size();


  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] ) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=4;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { x[j] = x[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { sigma[j] = sigma[j-a[2]];}
    for(j = a[3];j<n;j++) { lmd[j] = lmd[j-a[3]];}
  }



  int len;
  double pi;
  pi = 3.1415926535897932;
  len = n;
  NumericVector densitq(len),out(len);
  NumericVector I0(len),I(len), sign(len);


  for(int a=0;a<len;a++)
  {
    I0[a]=(x[a]<=mu[a]);
    I[a]= 1-I0[a];
    sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
    densitq[a] = sqrt(2/pi)*
      exp(-pow((x[a]-mu[a]),2)/(2*sigma[a]*sigma[a]*sign[a]))/
        ((1+lmd[a])*sigma[a]);
  }

  if(!logarithm)
  {
    for(int i = 0;i<len;i++)
    { out[i] = exp(densitq[i]);   }
  }
  else {out = densitq;}
  return out;
}
