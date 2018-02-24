#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector dsplitn(NumericVector x,NumericVector mu, NumericVector phi, NumericVector sigma, bool logarithm)
{
  int len;
  len = x.size();
  NumericVector density(len),densMean(len),densStd(len);
  NumericVector I0(len),I(len), sign(len);

    for(int a=0;a<len;a++)
    {
      I0[a]=(q[a]<=mu[a]);
      I[a]= 1-I0[a];
      sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
      if(q[a]<=mu[a])
      {
        densitylog[a] =2/(1+lmd[a])*R::pnorm5(q[a],mu[a],sigma[a],1,0);
      }
      else if(q[a]>mu[a])
      {
        densitylog[a] = 1/(1+lmd[a]) +
          2*lmd[a]/(1+lmd[a])*(R::pnorm5(q[a],mu[a],sigma[a],1,0)-1/2);
      }
    }

    if(!logarithm)
      {
        for(i = 0;i<n;i++){ out[i] = exp(densitylog[i]);   }
      }

    return out;

}


NumericVector psplitn(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int len;
  len = q.size();
  NumericVector density(len),densMean(len),densStd(len);
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
  return out;
}
