#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dmvtRcpp(NumericVector x, NumericVector sigma, bool log0, Function f) {
  NumericVector res = dmvt(x, sigma,type,df,log0);
  return res;
}
