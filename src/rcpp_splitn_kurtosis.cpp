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
NumericVector splitn_kurtosis(NumericVector lmd)
{
  int n,i,j;
  double pi;
  n =lmd.size();
  pi = 3.1415926535897932;

  NumericVector kurtosis(n);
  NumericVector k1(n),k2(n),k3;

  for(int i=0;i<n;i++){
    k1[i] = pow((lmd[i]-1),2);
    k2[i] = 8*(pi-3)*lmd[i]*lmd[i]+3*pow((pi-4),2)+8*(pi-3);
    k3[i] = pow((-2*pow((lmd[i]-1),2)+pi*lmd[i]*(lmd[i]-1)+pi),2);
    kurtosis[i] = k1[i]*k2[i]/k3[i];
  }
  return kurtosis;
}
