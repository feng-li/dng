#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector dsplitt(NumericVector x,NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd, bool logarithm)
{
  int a[5];
  int n,i,j;
  a[0] = x.size();
  a[1] = mu.size();
  a[2] = df.size();
  a[3] = phi.size();
  a[4] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] && a[0]==a[4] ) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=4;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { x[j] = x[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { df[j] = df[j-a[2]];}
    for(j = a[3];j<n;j++) { phi[j] = phi[j-a[3]];}
    for(j = a[4];j<n;j++) { lmd[j] = lmd[j-a[4]];}
  }

  NumericVector I0(n),I(n);
  NumericVector sign(n),densitylog(n),out(n);
  NumericVector lbeta0(n);

    for(i = 0;i<n;i++)
    {
      lbeta0[i] = ::Rf_lbeta(0.5*df[i],0.5);
      I0[i] = (x[i]<=mu[i]); // Logical values. 1, if y <= mu; 0, if y >mu.
      I[i] = (x[i]>mu[i]); //Logical values. 1, if y > mu; 0, if y <= mu.
      sign[i] = 1*I0[i]+lmd[i]*I[i]; // sign = 1 if y<=mu; sign = lmd.^2 if y>2
      densitylog[i] = (std::log(2)+(1+df[i])/2*(std::log(df[i])-std::log(df[i]+pow((-mu[i]+x[i]),2)/(pow(phi[i],2)*pow(sign[i],2))))-std::log(phi[i])-std::log(df[i])/2-lbeta0[i]-std::log(1+lmd[i]));
      out[i] = densitylog[i];
    }


  if(!logarithm)
  {
    for(i = 0;i<n;i++){ out[i] = exp(densitylog[i]);   }
  }

  return out;
}


