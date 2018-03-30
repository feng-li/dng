#include <Rcpp.h>
using namespace Rcpp;
//' Split-t distribution
//'
//' Density function for the split student-t distribution.
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
//' @aliases dng_dsplitt dsplitt
//' @param x vector of quantiles.
//' @param mu vector of location parameter. (The mode of the density)
//' @param df degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
//' @param phi vector of scale parameters (>0).
//' @param lmd vector of skewness parameters (>0). If is 1, reduced to
//' symmetric student t distribution.
//' @param logarithm logical; if TRUE, probabilities p are given as log(p).
//' @return \code{dsplitt} gives the density.  Invalid arguments will result in
//' return value NaN, with a warning.
//'
//' The numerical arguments other than n are recycled to the length of the
//' result. Only the first elements of the logical arguments are used.
//' @author Feng Li, Jiayue Zeng
//' @seealso \code{\link{splitt_mean}()},
//' \code{\link{splitt_var}()},\code{\link{splitt_skewness}()} and
//' \code{\link{splitt_kurtosis}()} for numerical characteristics of the
//' Split-t distribution.
//' @references Li, F., Villani, M., & Kohn, R. (2010). Flexible modeling of
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
//' dsplitt0 <- dsplitt(x, mu, df, phi, lmd, logarithm = TRUE)
//'
//' @export
// [[Rcpp::export]]
NumericVector dsplitt(NumericVector x,NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd, bool logarithm)
{
  int a[5];
  int n,i,j;
  a[0] = x.size();
  a[1] = mu.size();
  a[2] = df.size();
  a[3] = phi.size();
  a[4] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] && a[0]==a[4] ) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=4;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { x[j] = x[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { df[j] = df[j-a[2]];}
    for(j = a[3];j<n;j++) { phi[j] = phi[j-a[3]];}
    for(j = a[4];j<n;j++) { lmd[j] = lmd[j-a[4]];}
  }

  NumericVector I0(n),I(n);
  NumericVector sign(n),densitylog(n),out(n);
  NumericVector lbeta0(n);

    for(i = 0;i<n;i++)
    {
      lbeta0[i] = ::Rf_lbeta(0.5*df[i],0.5);
      I0[i] = (x[i]<=mu[i]); // Logical values. 1, if y <= mu; 0, if y >mu.
      I[i] = (x[i]>mu[i]); //Logical values. 1, if y > mu; 0, if y <= mu.
      sign[i] = 1*I0[i]+lmd[i]*I[i]; // sign = 1 if y<=mu; sign = lmd.^2 if y>2
      densitylog[i] = (std::log(2)+(1+df[i])/2*(std::log(df[i])-std::log(df[i]+pow((-mu[i]+x[i]),2)/(pow(phi[i],2)*pow(sign[i],2))))-std::log(phi[i])-std::log(df[i])/2-lbeta0[i]-std::log(1+lmd[i]));
      out[i] = densitylog[i];
    }


  if(!logarithm)
  {
    for(i = 0;i<n;i++){ out[i] = exp(densitylog[i]);   }
  }

  return out;
}
