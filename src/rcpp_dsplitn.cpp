#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dsplitn(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector lmd, bool logarithm)
{
  int a[4];
  int n,i,j;
  a[0] = x.size();
  a[1] = mu.size();
  a[2] = sigma.size();
  a[3] = lmd.size();


  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] ) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=4;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { x[j] = x[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { sigma[j] = sigma[j-a[2]];}
    for(j = a[3];j<n;j++) { lmd[j] = lmd[j-a[3]];}
  }



  int len;
  double pi;
  pi = 3.1415926535897932;
  len = n;
  NumericVector densitq(len),out(len);
  NumericVector I0(len),I(len), sign(len);


  for(int a=0;a<len;a++)
  {
    I0[a]=(x[a]<=mu[a]);
    I[a]= 1-I0[a];
    sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
    densitq[a] = sqrt(2/pi)*
      exp(-pow((x[a]-mu[a]),2)/(2*sigma[a]*sigma[a]*sign[a]))/
        ((1+lmd[a])*sigma[a]);
  }

  if(!logarithm)
  {
    for(int i = 0;i<len;i++)
    { out[i] = exp(densitq[i]);   }
  }
  else {out = densitq;}
  return out;
}



