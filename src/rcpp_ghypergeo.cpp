#include <Rcpp.h>
using namespace Rcpp;
NumericMatrix ghypergeo(NumericMatrix a, NumericMatrix b, NumericVector z,int k);
NumericMatrix pochhammer(NumericVector a, NumericVector n, LogicalVector log0);
// [[Rcpp::export]]
NumericMatrix ghypergeo(NumericMatrix a, NumericMatrix b, NumericVector z,int k) {
  int i,j,n,nObs,a_nCol,b_nCol, b_nRow;
  NumericVector n_series;
  k = 10;
  nObs   = a.nrow();
  a_nCol = a.ncol();
  b_nCol = b.ncol();
  b_nRow = b.nrow();
  n = k+1;
  for(i=0;i<n;i++)
  {  n_series = i;
  }
  NumericMatrix zpower_series(nObs, n);
  for(i=0;i<nObs;i++)
  {
    for(j=0;j<n;j++)
    {zpower_series(i,j)=n_series[j]*log(z[i]);
    }
  }

  NumericVector nfact_series;
  NumericMatrix ta(a_nCol,nObs), tb(b_nCol, b_nRow),an_series, bn_series;

  for(i=0;i<a_nCol;i++)
  {
    for(i=0; i<nObs; i++)
    {ta(i,j)=a(j,i); }
  }

  for(i=0;i<b_nCol;i++)
  {
    for(i=0; i<b_nRow; i++)
    {tb(i,j)=b(j,i); }
  }


  for(i=0;i<n;i++)
  {
    nfact_series[i] = R::lgammafn(n_series[i]+1) ;
  }
  //an_series    = pochhammer(ta, n_series, TRUE);
  //bn_series    = pochhammer(tb, n_series, TRUE);

  NumericMatrix out(n,1);
  for(i=0;i<n;i++)
  {
    out(i,0)=i;
  }
  return (out);

}


