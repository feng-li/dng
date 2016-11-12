#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector vech2mRcpp(NumericVector vech, bool diag, Function vech2m) {
  NumericVector res = vech2m(vech, diag);
  return res;
}
