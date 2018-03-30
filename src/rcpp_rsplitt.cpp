#include <Rcpp.h>
using namespace Rcpp;
//' random generation for Split-t distribution
//'
//' Density, distribution function, quantile function and random generation for
//' the split student-t distribution.
//'
//' The random variable y follows a split-t distribution with \eqn{\nu}>0
//' degrees of freedom, y~t(\eqn{\mu}, \eqn{\phi}, \eqn{\lambda}, \eqn{\nu}),
//' if its density function is of the form
//'
//' \deqn{C K(\mu, \phi, \nu,)I(y\leq\mu) + C K(\mu, \lambda \phi,
//' \nu)I(y>\mu), } where, \deqn{K(\mu, \phi, \nu,) =[\nu/(\nu+(y-\mu)^2 /\phi
//' ^2)]^{(\nu+1)/2} } is the kernel of a student \eqn{t} density with variance
//' \eqn{\phi ^2\nu/(\nu-2)} and \deqn{c = 2[(1+\lambda)\phi (\sqrt \nu)
//' Beta(\nu/2,1/2)]^{-1} }is the normalization constant.
//'
//' @aliases rsplitt dng_rsplitt
//' @param n number of observations. If length(n) > 1, the length is taken to
//' be the number required.
//' @param mu vector of location parameter. (The mode of the density)
//' @param df degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
//' @param phi vector of scale parameters (>0).
//' @param lmd vector of skewness parameters (>0). If is 1, reduced to
//' symmetric student t distribution.
//' @return \code{rsplitt} generates random deviates. Invalid arguments will
//' result in return value NaN, with a warning.
//'
//' The length of the result is determined by \code{n} for \code{rsplitt}.
//'
//' The numerical arguments other than n are recycled to the length of the
//' result. Only the first elements of the logical arguments are used.
//' @author Feng Li, Jiayue Zeng
//' @seealso \code{\link{splitt_mean}()},
//' \code{\link{splitt_var}()},\code{\link{splitt_skewness}()} and
//' \code{\link{splitt_kurtosis}()} for numerical characteristics of the
//' Split-t distribution.
//' @references Li, F., Villani, M., & Kohn, R. (2009). Flexible modeling of
//' conditional distributions using smooth mixtures of asymmetric student t
//' densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
//' @keywords distribution asymmetric student-t
//' @examples
//'
//' n <- 5
//' x <- c(0.25,0.5,0.75)
//' q <- c(0.25,0.5,0.75)
//' p <- c(0.25,0.5,0.75)
//' mu <- c(0,1,2)
//' df <- rep(10,3)
//' phi <- c(0.5,1,2)
//' lmd <- c(1,2,3)
//'
//' rsplitt0 <- rsplitt(n, mu, df, phi, lmd)
//' @export
// [[Rcpp::export]]
NumericVector rsplitt(int n, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  int i,j;
  NumericVector u(n),out(n);
  for(i = 0; i<n; i++)
  {
    u[i] = R::runif(0,1);
  }

 //out = qsplitt(u, mu, df, phi, lmd);
  int a[5];
  a[0] = u.size();
  a[1] = mu.size();
  a[2] = df.size();
  a[3] = phi.size();
  a[4] = lmd.size();

  for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
  for(j = a[2];j<n;j++) { df[j] = df[j-a[2]];}
  for(j = a[3];j<n;j++) { phi[j] = phi[j-a[3]];}
  for(j = a[4];j<n;j++) { lmd[j] = lmd[j-a[4]];}

  NumericVector mu_long(n),df_long(n),phi_long(n),lmd_long(n);
  NumericVector I0(n),I(n);
  NumericVector p0std(n), y0std(n);



  for(i = 0;i<n;i++)
  {
    I0[i] = (u[i]<=(1/(1+lmd[i])));

    if(I0[i])
    {
      p0std[i] = u[i]*(1+lmd[i])/2;
      y0std[i] = R::qt(p0std[i], df[i], TRUE, FALSE);
      out[i] = y0std[i]*phi[i]+mu[i] ;
    }

    else
    {
      p0std[i] = (u[i]-1/(1+lmd[i]))*(1+lmd[i])/(2*lmd[i])+0.5;
      y0std[i] = R::qt(p0std[i], df[i], TRUE, FALSE);
      out[i] = y0std[i]*(phi[i]*lmd[i])+mu[i] ;
    }

  }




  return out;
}
