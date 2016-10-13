#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector logDensFun(NumericMatrix u, NumericVector theta , NumericVector delta)
{
  int i,j,n,m;
  n = u.nrow();
  m = u.ncol();

  NumericVector L5(n), L6(n), logCplDensObs(n), out_log(n), rs1(n), rs2(n);
  NumericMatrix TC1(n,m), TC2_log(n,m);

// The density function
  for(i= 0;i<n;i++)
  {
    for(j=0;j<m;j++)
    {
      TC1(i,j) = (1-pow((1-u(i,j)),theta[i])) ;// FIXME: Numerically instable if theta -> Inf,  then TC1-> 1
      TC2_log(i,j) = (-1+theta[i])*log(1-u(i,j)); // TC2 = (1-u)^(-1+theta)
    }
  }

for(i=0;i<n;i++)
{
  L5[i]  = -1;
  rs1[i] = 0;
  rs2[i] = 0;
  for(j=0;j<m;j++)
  {
    L5[i]  = L5[i]  + pow(TC1(i,j),(-delta[i]));
    rs1[i] = rs1[i] + log(TC1(i,j));
    rs2[i] = rs2[i] + TC2_log(i,j);
  }
}

for(i=0;i<n;i++)
{
  L6[i] = 1-pow(L5[i],(-1/delta[i])); // FIXME: log(L6)->Inf when u->1,  v->1
}

for(i=0;i<n;i++)
{
  logCplDensObs[i] = ((-1-delta[i])*rs1[i]+rs2[i]  -
    2*(1+delta[i])/delta[i]*log(L5[i])+
    (-2+1/theta[i])*log(L6[i])+
    log(-1+theta[i]+pow(L5[i],(1/delta[i]))*L6[i]*(1+delta[i])*theta[i]));
}

    //out_log = matrix(logCplDensObs);

    return logCplDensObs;

}
