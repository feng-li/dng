#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix FUN(int i, NumericMatrix x, NumericVector mu, NumericVector df, NumericMatrix rho, NumericVector uIdx )
{
  Environment Callfunc("package:grand");
  Function mult = Callfunc["mult"];
  Function vech2m = Callfunc["vech2m"];
  Environment base("package:base");
  Function solve = base["solve"];
  Function t = base["t"];

  NumericMatrix out;
  int rho_ncol = rho.ncol();
  int uIdx_n = uIdx.size();
  NumericVector rhoi(rho_ncol);
  for(int j=0;j<rho_ncol;j++)
  {    rhoi[j] = rho(i,j);  }
  NumericMatrix Sigma0 = vech2m(rhoi, 0);
  NumericMatrix Sigma;
  for(int a=0;a<uIdx_n;a++)
  {
    for(int b=0;b<uIdx_n;b++)
    {
      Sigma(a,b) = Sigma0(uIdx[a],uIdx[b]);
      if(TRUE ==!R_FINITE(Sigma(a,b)) | Sigma(a,b)<0)
      { out = 0; }
    }
  }

  int p = Sigma.nrow();
  int q = Sigma.ncol();
  double v = df[i];
  NumericVector x_i;
  int x_ncol = x.ncol();
  for(int j=0;j<x_ncol;j++)
  {
    x_i[j] = x(i,j);
  }

  NumericVector C2 = mult(mult(t(x_i-mu),solve(Sigma)),(x_i-mu));
  NumericVector solveSximu = mult(solve(Sigma),(x_i-mu));
  for(int a;a<p;a++)
  {
    for(int b;b<q;b++)
    {
      out(a,b) = -(v+p)/2 * pow((1+C2[b]/v),(-1)) * 1/v*(2*solveSximu(a,b));
    }
  }

  return out;

}
