#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dmvNormVecFun(int i,NumericVector x, NumericVector rho)
{
  Environment callmvtnorm("package:mvtnorm");
  Function dmvnorm = callmvtnorm["mvtnorm"];
  Environment Callfunc("package:dng");
  Function vech2m = Callfunc["vech2m"];
  NumericMatrix xi;
  int n=sizeof(rho);
  for(int j=0;j<n;j++)
  {
    xi(1,j) = x(i,j);
    }
  NumericVector Sigma = vech2m(rho, 0);
  NumericVector out = dmvnorm(xi, Sigma, 1);

  return out;
}
