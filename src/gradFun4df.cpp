#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector gradFun4df(int i, NumericMatrix rho, NumericVector df , NumericMatrix u_quantile)
{
  Environment Callfunc("package:grand");
  Function mult = Callfunc["mult"];
  Function vech2m = Callfunc["vech2m"];
  Environment base("package:base");
  Function solve = base["solve"];
  Function t = base["t"];

  int j;
  int rho_ncol = rho.ncol();
  //int u_quantile_ncol = u_quantile.ncol();
  NumericVector rhoi(rho_ncol);
  for(j=0;j<rho_ncol;j++)
  {    rhoi[j] = rho(i,j);  }
  NumericMatrix Sigma(1,rho_ncol);
  Sigma = vech2m(rhoi, 0);

  double v = df[i];


  NumericMatrix x(1,u_quantile_ncol);
  for(j=0;j<u_quantile_ncol;j++)
  { x(0,j) = u_quantile(i,j);  }
  int mu = 0;
  int q = Sigma.nrow();

  //C2 = as.vector(t(x-mu)%*%solve(Sigma)%*%(x-mu))
  NumericVector C2 = mult(t(x-mu),solve(Sigma, (x-mu)));

  NumericVector out= ((C2-q -(C2+v)*log((C2+v)/v)+
  (C2+v)*(-R::digamma(v/2)+R::digamma((q+v)/2)))/(2*(C2+v)));


  return out;

}
