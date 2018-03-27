#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector splitn_mean(NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int a[3];
  int n,i,j;
  double pi;
  pi = 3.1415926535897932;
  a[0] = mu.size();
  a[1] = sigma.size();
  a[2] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2]) {n = a[0];}
  else
  {
    n=a[0];
    for(i = 1;i<=2;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { mu[j] = mu[j-a[0]];}
    for(j = a[1];j<n;j++) { sigma[j] = sigma[j-a[1]];}
    for(j = a[2];j<n;j++) { lmd[j] = lmd[j-a[2]];}
  }

  NumericVector mean(n);

  for(int i=0;i<n;i++){
    mean[i] = sqrt(2/pi)*(lmd[i]-1)*sigma[i]+mu[i];
  }
  return mean;
}
