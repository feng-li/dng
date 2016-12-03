#include <Rcpp.h>
using namespace Rcpp;
/*** R
vech2m <- function(vech, diag = TRUE)
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
NumericVector gradFun4df(int i, NumericMatrix rho, NumericVector df , NumericVector u_quantile)
{
  Rcpp::Environment G = Rcpp::Environment::global_env();
  Rcpp::Function vech2m = G["vech2m"];
  int rho_nrow, rho_ncol,j;
  rho_nrow = rho.nrow();
  rho_ncol = rho.ncol();
  NumericVector rho_i(rho_ncol),u_quantile_i(rho_ncol);
  for(j=0;j<rho_ncol;j++)
  {
    rho_i[j] = rho(i,j);
    u_quantile_i[j] = u_quantile(i,j);
  }
  NumericMatrix Sigma = vech2m(rho_i, FALSE);
  double v = df[i];
  NumericMatrix x = Rcpp::as<NumericMatrix>(u_quantile_i); // col-vector
  int mu = 0;
  int q = Sigma.nrow();
//C2 = as.vector(t(x-mu)%*%solve(Sigma)%*%(x-mu))

  int n;
  NumericVector C2(rho_ncol), out(rho_ncol);
  NumericMatrix solvex = R::solve(Sigma,(x-mu));

  for(n=0;n<q;n++)
  {
    C2[n] = 0;
    for(j=0;j<rho_ncol;j++)
    {
      C2[n] = C2[n]+(x[j]-mu)*solvex(j,n);
    }
  }

  for(j=0;j<rho_ncol;j++)
    {
      out[j] = ((C2[j]-q -(C2[j]+v)*log((C2[j]+v)/v)+
      (C2[j]+v)*(-R::digamma(v/2)+R::digamma((q+v)/2)))/(2*(C2[j]+v)));
    }

    return out;

}
