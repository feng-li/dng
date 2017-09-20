#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rsplitt(int n, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  int i,j;
  NumericVector u(n),out(n);
  for(i = 0; i<n; i++)
  {
    u[i] = R::runif(0,1);
  }

 //out = qsplitt(u, mu, df, phi, lmd);
  int a[5];
  a[0] = u.size();
  a[1] = mu.size();
  a[2] = df.size();
  a[3] = phi.size();
  a[4] = lmd.size();

  for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
  for(j = a[2];j<n;j++) { df[j] = df[j-a[2]];}
  for(j = a[3];j<n;j++) { phi[j] = phi[j-a[3]];}
  for(j = a[4];j<n;j++) { lmd[j] = lmd[j-a[4]];}

  NumericVector mu_long(n),df_long(n),phi_long(n),lmd_long(n);
  NumericVector I0(n),I(n);
  NumericVector p0std(n), y0std(n);



  for(i = 0;i<n;i++)
  {
    I0[i] = (u[i]<=(1/(1+lmd[i])));

    if(I0[i])
    {
      p0std[i] = u[i]*(1+lmd[i])/2;
      y0std[i] = R::qt(p0std[i], df[i], TRUE, FALSE);
      out[i] = y0std[i]*phi[i]+mu[i] ;
    }

    else
    {
      p0std[i] = (u[i]-1/(1+lmd[i]))*(1+lmd[i])/(2*lmd[i])+0.5;
      y0std[i] = R::qt(p0std[i], df[i], TRUE, FALSE);
      out[i] = y0std[i]*(phi[i]*lmd[i])+mu[i] ;
    }

  }




  return out;
}



