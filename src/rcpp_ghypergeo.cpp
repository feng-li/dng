#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ghypergeo(NumericMatrix a, NumericMatrix b, NumericVector z,int k) {
  int i,j,n,m,nObs,a_nCol,b_nCol, b_nRow;
  NumericVector n_series;
  k = 10;
  nObs   = a.nrow();
  a_nCol = a.ncol();
  b_nCol = b.ncol();
  b_nRow = b.nrow();
  n = k+1;

  NumericVector out(nObs) ;
  for(i=0;i<nObs;i++) { out[i]=0; }

  for(i=0;i<n;i++)
  {  n_series = i;
  }
  NumericMatrix zpower_series(nObs, n),tzpower_series(n, nObs), nlog(n, nObs);
  for(i=0;i<nObs;i++)
  {
    for(j=0;j<n;j++)
    {zpower_series(i,j) = n_series[j]*log(z[i]);
    }
  }


  for(i=0;i<n;i++)
  {
    for(j=0;j<nObs;j++)
    {tzpower_series(i,j) = zpower_series(j,i);
    }
  }



  NumericVector nfact_series;
  NumericMatrix ta(a_nCol,nObs), tb(b_nCol, b_nRow),an_series(a_nCol,n), bn_series(b_nCol,n);

  for(i=0;i<a_nCol;i++)
  {
    for(j=0; j<nObs; j++)
    {ta(i,j)=a(j,i); }
  }

  for(i=0;i<b_nCol;i++)
  {
    for(j=0; j<b_nRow; j++)
    {tb(i,j)=b(j,i); }
  }


  for(i=0;i<n;i++)
  {
    nfact_series[i] = R::lgammafn(n_series[i]+1) ;
  }

  for(i=0;i<a_nCol;i++)
  {
    for(j=0;j<n;j++)
    {
      an_series(i,j) = R::lgammafn(ta(i,j)+n_series[i]) - R::lgammafn(ta(i,j));
    }
  }

  for(i=0;i<b_nCol;i++)
  {
    for(j=0;j<n;j++)
    {
      bn_series(i,j) = R::lgammafn(tb(i,j)+n_series[i]) - R::lgammafn(tb(i,j));
    }
  }

  NumericMatrix an(n,nObs),bn(n,nObs),tan_series(n,a_nCol),tbn_series(n,b_nCol);

  for(i=0;i<n;i++)
  {
    for(j=0; j<a_nCol; j++)
    {tan_series(i,j)=an_series(j,i); }
  }

  for(i=0;i<n;i++)
  {
    for(j=0; j<b_nCol; j++)
    {tbn_series(i,j) = bn_series(j,i); }
  }


  an=0;
  bn=0;
  for(i=0;i<n;i++)
  {
    for(j=0;j<nObs;j++)
    {
      for(m=0;m<a_nCol;m++)
      {
        an(i,j) = an(i,j)+tan_series(i,m);
      }
    }
  }

  for(i=0;i<n;i++)
  {
    for(j=0;j<nObs;j++)
    {
      for(m=0;m<b_nCol;m++)
      {
        bn(i,j) = bn(i,j)+tbn_series(i,m);
      }
    }
  }



  for(i=0;i<n;i++)
  {
    for(j=0;j<nObs;j++)
    {
      nlog(i,j) = an(i,j) - bn(i,j) + tzpower_series(i,j) - nfact_series[i];
    }
  }

  for(i=0;i<nObs;i++)
  {
    for(j=0;j<n;j++)
    {
      out[i] = out[i]+exp(nlog(j,i));
    }
  }


  return (out);

}


