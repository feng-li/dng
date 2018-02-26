#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector psplitn(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int len;
  len = q.size();
  NumericVector densitq(len),out(len);
  NumericVector I0(len),I(len), sign(len);


  for(int a=0;a<len;a++)
  {
    I0[a]=(q[a]<=mu[a]);
    I[a]= 1-I0[a];
    sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
    densitq[a] = sqrt(2/3.1415926)*
      exp(-pow((q[a]-mu[a]),2)/(2*sigma[a]*sigma[a]*sign[a]))/
        ((1+lmd[a])*sigma[a]);
  }

  out = densitq;
  return out;
}
