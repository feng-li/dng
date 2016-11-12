#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dmvNormVecFunRcpp(int i,NumericVector x, NumericVector rho, Function f) {
  NumericVector res = dmvNormVecFun(i, x, rho);
  return res;
}
