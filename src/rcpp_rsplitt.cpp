#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rsplitt(int n, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
extern NumericVector qsplitt(NumericVector p,NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);

NumericVector rsplitt(int n, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  NumericVector u(n),out(n);
  for(int i = 0; i<n; i++)
  {
    u[i] = R::runif(0,1);//!!!
  }

  out = qsplitt(u, mu, df, phi, lmd);


  return out;
}



