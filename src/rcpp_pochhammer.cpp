#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pochhammer(NumericVector a, IntegerVector n, LogicalVector log0) {
  int i,j,a_len, n_len;
  a_len = a.size();
  n_len = n.size();
  NumericMatrix out_log(a_len,n_len),out(a_len,n_len);
  NumericMatrix a_mat(a_len,n_len),n_mat(a_len,n_len);
  for(i = 0;i<a_len;i++)
  {
    for(j = 0;j<n_len;j++)
    {
      a_mat(i,j) = a[i];
      n_mat(i,j) = n[j];
    }
  }

  for(i=0;i<a_len;i++)
  {
    for(j=0;j<n_len;j++)
    {out_log(i,j) = R::lgammafn(a_mat(i,j)+n_mat(i,j))-R::lgammafn(a_mat(i,j));}
  }


  if(log0[0])
  {
    out = out_log;
  }
  else
  {
    for(i=0;i<a_len;i++)
    {
      for(j=0;j<n_len;j++)
      {
        out(i,j) = exp(out_log(i,j));
      }
    }

  }
  return(out);

}


