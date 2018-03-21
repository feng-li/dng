#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector qsplitn(NumericVector p,NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int a[4];
  int n,i,j;
  a[0] = p.size();
  a[1] = mu.size();
  a[2] = sigma.size();
  a[3] = lmd.size();


  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] ) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=3;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { p[j] = p[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { sigma[j] = sigma[j-a[2]];}
    for(j = a[3];j<n;j++) { lmd[j] = lmd[j-a[3]];}
  }

  NumericVector p0(n),quantile(n);

  for(int i=0;i<n;i++)
  {
    if(p[i]<=(1/(1+lmd[i])))
    {
      p0[i] = (1+lmd[i])*p[i]/2;
      quantile[i] = R::pnorm5(p0[i],mu[i],sigma[i],1,0);

    }

    else
    {
      p0[i] = (p[i]-(1-lmd[i])/(1+lmd[i]))*(1+lmd[i])/(2*lmd[i]);
      quantile[i] = R::pnorm5(p0[i],mu[i],(sigma[i]*lmd[i]),1,0);
    }

  }

  return quantile;

}
