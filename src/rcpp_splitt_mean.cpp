#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector splitt_mean(NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  int a[4];
  int n,i,j;
  a[0] = mu.size();
  a[1] = df.size();
  a[2] = phi.size();
  a[3] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3]) {n = a[0];}
  else
  {
    n=a[0];
    for(i = 1;i<=3;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { mu[j] = mu[j-a[0]];}
    for(j = a[1];j<n;j++) { df[j] = df[j-a[1]];}
    for(j = a[2];j<n;j++) { phi[j] = phi[j-a[2]];}
    for(j = a[3];j<n;j++) { lmd[j] = lmd[j-a[3]];}
  }

  NumericVector h(n),mean(n);
  NumericVector beta0(n);

  for(int i=0;i<n;i++){
    beta0[i]=R::beta(df[i]*0.5,0.5);
    h[i] = 2*pow(df[i],0.5)*phi[i]*(lmd[i]-1)/((df[i]-1)*beta0[i]);
    mean[i] = mu[i]+h[i];
  }
  return mean;
}
