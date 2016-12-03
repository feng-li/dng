#include <Rcpp.h>
using namespace Rcpp;
/*** R
dmvNormVecFunR <- function(i, x, rho)
{
  require("mvtnorm")
  Sigma = vech2m(rho[i, ], diag = FALSE)
  out <- dmvnorm(x = x[i, , drop = FALSE],
                 sigma = Sigma,
                 log = TRUE)
  return(out)
}
*/

// [[Rcpp::export]]
NumericVector dmvNormVecFun(int i,NumericVector x, NumericVector rho) {

  Rcpp::NumericVector res;
  Rcpp::Environment G = Rcpp::Environment::global_env();
  Rcpp::Function inte = G["dmvNormVecFunR"];
  res= dmvNormVecFun(i, x, rho);
  return res;
}
