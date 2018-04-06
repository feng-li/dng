#include <Rcpp.h>
using namespace Rcpp;
//' Distribution function of Split-t distribution
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
//' @aliases psplitt dng_psplitt
//' @param q vector of quantiles.
//' @param mu vector of location parameter. (The mode of the density)
//' @param df degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
//' @param phi vector of scale parameters (>0).
//' @param lmd vector of skewness parameters (>0). If is 1, reduced to
//' symmetric student t distribution.
//' @return \code{psplitt} gives the distribution function.  (\code{dsplitt},
//' \code{psplitt}, \code{qsplitt} and \code{rsplitt} are all vectors.)
//'
//' Invalid arguments will result in return value NaN, with a warning.
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
//' psplitt0 <- psplitt(q, mu, df, phi, lmd)
//'
//' @export
// [[Rcpp::export]]
NumericVector psplitt(NumericVector q, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
extern double ibeta(double x, double a, double b, bool log0, bool reg);

NumericVector psplitt(NumericVector q, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  double ibeta0;
  int a[5];
  int n,i,j;
  a[0] = q.size();
  a[1] = mu.size();
  a[2] = df.size();
  a[3] = phi.size();
  a[4] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] && a[0]==a[4]) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=5;i++)   { if(a[i]>n) n = a[i];}

    for(j = a[0];j<n;j++) { q[j] = q[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { df[j] = df[j-a[2]];}
    for(j = a[3];j<n;j++) { phi[j] = phi[j-a[3]];}
    for(j = a[4];j<n;j++) { lmd[j] = lmd[j-a[4]];}
  }

  NumericVector I0(n),I(n), sign(n), sign2(n);
  NumericVector A(n),BetaRegUpper(n);
  NumericVector out(n);

  for(i = 0;i<n;i++){
    I0[i] = (q[i]<=mu[i]);
    I[i]  = 1-I0[i];
    sign[i]  = 1*I0[i]+lmd[i]*I[i];
    sign2[i] = -1*I0[i]+1*I[i];

    A[i] = df[i]*pow(sign[i],2)*pow(phi[i],2)/(df[i]*pow(sign[i],2)*pow(phi[i],2)+pow((q[i]-mu[i]),2));

    ibeta0 = ibeta(A[i], df[i]*0.5, 0.5, FALSE, TRUE);
    BetaRegUpper[i] = 1-ibeta0;

    out[i] = (1/(1+lmd[i]) + sign[i]*sign2[i]/(1+lmd[i])*BetaRegUpper[i]);

  }

  return out;
}
