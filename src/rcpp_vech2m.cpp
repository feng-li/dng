#include <Rcpp.h>
using namespace Rcpp;
/*** R
vech2m <- function(vech, diag)
{
  c <- length(vech)
  p <- (1 - 2*diag +sqrt(8*c+1))/2 # n^2 -n  =  2c; or n^2 + n  =  2c
  if(!identical(p,  round(p)))
  {
    stop ("Input vector length not match with the length of matrix triangular.")
  }
  out <- matrix(1, p, p)
    idx.lower <- as.vector(lower.tri(out, diag = diag))
    out[idx.lower] <- vech
    idx.upper <- as.vector(upper.tri(out, diag = diag))
    out.t <- t(out)
    out[idx.upper] <- out.t[idx.upper]
  return(out)
}
*/

// [[Rcpp::export]]
double vech2m(Rcpp::NumericVector vech, bool diag)
{
  Rcpp::NumericVector res;
  Rcpp::Environment G = Rcpp::Environment::global_env();
  Rcpp::Function vech2m = G["vech2m"];
  res = vech2m(vech, diag);
  return res[0];
}
