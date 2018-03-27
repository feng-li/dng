#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector psplitn(NumericVector q,NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int a[4];
  int n,i,j;
  a[0] = q.size();
  a[1] = mu.size();
  a[2] = sigma.size();
  a[3] = lmd.size();


  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] ) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=4;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { q[j] = q[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { sigma[j] = sigma[j-a[2]];}
    for(j = a[3];j<n;j++) { lmd[j] = lmd[j-a[3]];}
  }

  int len;
  len = n;
  NumericVector density(len), out(len);
  NumericVector I0(len),I(len), sign(len);

  for(int a=0;a<len;a++)
  {
    I0[a]=(q[a]<=mu[a]);
    I[a]= 1-I0[a];
    sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
    if(q[a]<=mu[a])
    {
      out[a] =2/(1+lmd[a])*R::pnorm5(q[a],mu[a],sigma[a],1,0);
    }
    else if(q[a]>mu[a])
    {
      out[a] = (1-lmd[a])/(1+lmd[a]) +
        2*lmd[a]/(1+lmd[a])*(R::pnorm5(q[a],mu[a],sigma[a]*lmd[a],1,0)-1/2);
    }
  }

  return out;

}
