#include <Rcpp.h>
using namespace Rcpp;
double dmvNormVecFun(int i,NumericVector x, NumericVector rho);
double vech2m(Rcpp::NumericVector vech, bool diag);

// [[Rcpp::export]]
double dmvNormVecFun(int i,NumericVector x, NumericVector rho)
{

  Environment callmvtnorm("package:mvtnorm");
  Function dmvnorm = callmvtnorm["mvtnorm"];

  NumericVector Sigma = vech2m(rho, 0);
  double out = dmvnorm(x, Sigma, 1);

  return out;
}
