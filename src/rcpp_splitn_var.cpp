#include <Rcpp.h>
using namespace Rcpp;
//' Moments of the Split Normal distribution
//'
//' Compute the mean, variance, skewness and kurtosis for the split normal
//' distribution.
//'
//' The random variable y follows a split-normal distribution, which has
//' density: \deqn{1/(1+\lambda)\sigma \sqrt(2/\pi) exp{-(y-\mu)*2/2\sigma^2},
//' if y<=\mu} \deqn{1/(1+\lambda)\sigma \sqrt(2/\pi) exp{-(y-\mu)*2/2\sigma^2
//' \lambda^2}, if y>\mu} where \eqn{\sigma>0} and \eqn{\lambda>0}. The
//' Split-normal distribution reduce to normal distribution when \eqn{\lmd=1}.
//'
//' @aliases splitn_mean splitn_var splitn_skewness splitn_kurtosis
//' dng_splitn_mean dng_splitn_var dng_splitn_skewness dng_splitn_kurtosis
//' @param mu vector of location parameter. (The mode of the density)
//' @param sigma vector of standard deviations.
//' @param lmd vector of skewness parameters (>0). If is 1, reduce to normal
//' distribution.
//' @return \code{splitn_mean} gives the mean.  \code{splitn_var} gives the
//' variance.  \code{splitn_skewness} gives the skewness.
//' \code{splitn_kurtosis} gives the kurtosis.  (\code{splitn_mean},
//' \code{splitn_var},\code{splitn_skeness} and \code{splitn_kurtosis} are all
//' vectors.
//' @author Feng Li, Jiayue Zeng
//' @seealso \code{\link{psplitn}()} \code{\link{dsplitn}()} and for the
//' split-normal distribution.
//' @references Li, F., Villani, M., & Kohn, R. (2010). Flexible modeling of
//' conditional distributions using smooth mixtures of asymmetric student t
//' densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
//' and Villani, M., & Larsson, R. (2006) The Multivariate Split Normal
//' Distribution and Asymmetric Principal Components Analysis. Sveriges
//' Riksbank Working Paper Series, No. 175.
//' @keywords distribution asymmetric normal
//' @examples
//'
//' mu <- c(0,1,2)
//' sigma <- c(0.5,1,2)
//' lmd <- c(1,2,3)
//'
//' mean0 <- splitn_mean(mu, sigma, lmd)
//' var0 <- splitn_var(sigma, lmd)
//' skewness0 <- splitn_skewness(sigma, lmd)
//' kurtosis0 <- splitn_kurtosis(lmd)
//' @export
// [[Rcpp::export]]
NumericVector splitn_var(NumericVector sigma, NumericVector lmd)
{
  int a[2];
  int n,i,j;
  double pi;
  pi = 3.1415926535897932;
  a[0] = sigma.size();
  a[1] = lmd.size();

  if(a[0]==a[1]) {n = a[0];}
  else
  {
    n=a[0];
    if(a[1]>n) n = a[1];
    for(j = a[0];j<n;j++) { sigma[j] = sigma[j-a[0]];}
    for(j = a[1];j<n;j++) { lmd[j] = lmd[j-a[1]];}
  }

  NumericVector var(n),b(n);

  for(int i=0;i<n;i++){
    b[i] = -2/pi*pow((lmd[i]-1),2)+lmd[i]*(lmd[i]-1)+1;
    var[i] = b[i]*pow(sigma[i],2);
  }
  return var;
}
