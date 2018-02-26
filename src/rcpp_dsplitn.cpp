#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector dsplitn(NumericVector x,NumericVector mu, NumericVector phi, NumericVector sigma, bool logarithm)
{
  int len;
  len = x.size();
  NumericVector densitylog(len), out(len);
  NumericVector I0(len),I(len), sign(len);

    for(int a=0;a<len;a++)
    {
      I0[a]=(x[a]<=mu[a]);
      I[a]= 1-I0[a];
      sign[a]=1*I0[a]+phi[a]*phi[a]*I[a];
      if(x[a]<=mu[a])
      {
        densitylog[a] =2/(1+phi[a])*R::pnorm5(x[a],mu[a],sigma[a],1,0);
      }
      else if(x[a]>mu[a])
      {
        densitylog[a] = 1/(1+phi[a]) +
          2*phi[a]/(1+phi[a])*(R::pnorm5(x[a],mu[a],sigma[a],1,0)-1/2);
      }
    }

    if(!logarithm)
      {
        for(int i = 0;i<len;i++)
          { out[i] = exp(densitylog[i]);   }
      }

    return out;

}
