#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitt Percentile for the split-t distribution.
//' @export
// [[Rcpp::export]]
NumericVector psplitt(NumericVector q, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  double ibeta0;
  double pbeta0;
  int a[5];
  int n,i,j;
  a[0] = q.size();
  a[1] = mu.size();
  a[2] = df.size();
  a[3] = phi.size();
  a[4] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3] && a[0]==a[4]) {n = a[0];}
  else
  {
    n = a[0];
    for(i = 1;i<=4;i++)   { if(a[i]>n) n = a[i];}

    for(j = a[0];j<n;j++) { q[j] = q[j-a[0]];}
    for(j = a[1];j<n;j++) { mu[j] = mu[j-a[1]];}
    for(j = a[2];j<n;j++) { df[j] = df[j-a[2]];}
    for(j = a[3];j<n;j++) { phi[j] = phi[j-a[3]];}
    for(j = a[4];j<n;j++) { lmd[j] = lmd[j-a[4]];}
  }

  NumericVector I0(n),I(n), sign(n), sign2(n);
  NumericVector A(n),BetaRegUpper(n);
  NumericVector out(n);

  for(i = 0;i<n;i++){
    I0[i] = (q[i]<=mu[i]);
    I[i]  = 1-I0[i];
    sign[i]  = 1*I0[i]+lmd[i]*I[i];
    sign2[i] = -1*I0[i]+1*I[i];

    A[i] = df[i]*pow(sign[i],2)*pow(phi[i],2)/(df[i]*pow(sign[i],2)*pow(phi[i],2)+pow((q[i]-mu[i]),2));

    pbeta0 = ::Rf_pbeta(A[i],df[i]*0.5, 0.5, TRUE, TRUE);
    ibeta0 = exp(pbeta0);

    BetaRegUpper[i] = 1-ibeta0;

    out[i] = (1/(1+lmd[i]) + sign[i]*sign2[i]/(1+lmd[i])*BetaRegUpper[i]);

  }

  return out;
}
