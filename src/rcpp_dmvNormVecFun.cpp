#include <Rcpp.h>
using namespace Rcpp;
NumericVector dmvNormVecFun(int i,NumericVector x, NumericVector rho);
double vech2m(Rcpp::NumericVector vech, bool diag);
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
NumericVector dmvNormVecFun(int i,NumericVector x, NumericMatrix rho)
{

  Environment callmvtnorm("package:mvtnorm");
  Function dmvnorm = callmvtnorm["mvtnorm"];
  NumericVector rhoi;
  int rho_ncol = rho.ncol();
  for(int j=0;j<rho_ncol;j++)
  {    rhoi[j] = rho(i,j);  }
  NumericVector Sigma = vech2m(rhoi, FALSE);
  NumericVector out = dmvnorm(x = x[i, , drop = FALSE],
                               sigma = Sigma,
                               log = TRUE)

  return out;
}
