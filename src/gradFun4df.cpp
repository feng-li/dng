#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gradFun4df(int i, NumericMatrix rho, NumericVector df , NumericVector u_quantile)
{
  Sigma <- vech2m(rho[i, ], diag = FALSE)
  v <- df[i]
  x <- matrix(u.quantile[i, ]) # col-vector
  mu <- 0
  q <- dim(Sigma)[1]
//C2 <- as.vector(t(x-mu)%*%solve(Sigma)%*%(x-mu))
  C2 <- as.vector(t(x-mu)%*%solve(Sigma, (x-mu)))

    out <- ((C2-q -(C2+v)*log((C2+v)/v)+
      (C2+v)*(-digamma(v/2)+digamma((q+v)/2)))/(2*(C2+v)))


  return out;

}
