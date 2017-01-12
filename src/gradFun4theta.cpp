#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gradFun4theta(NumericMatrix u, NumericVector theta , NumericVector delta)
{
  int i,j,n,m;
  n = u.nrow();
  m = u.ncol();

  NumericVector outtheta(n),L1(n), SD12(n), SD34(n), C1(n), C2(n), rowsum1(n), rowsum2(n);
  NumericMatrix T1(n,m), PT1(n,1);
  double la;



  for(i= 0;i<n;i++)
  {
    L1[i] = -1;
    for(j= 0;j<m;j++)
    {
      T1(i,j)  = 1 - pow((1-u(i,j)),theta[i]);
      L1[i]    = L1[i]   + pow(T1(i,j),(-delta[i]));
      SD12[i]  = SD12[i] + pow(T1(i,j),(1+delta[i]))*pow((1-u(i,j)),theta[i])*log(1-u(i,j));
      SD34[i]  = SD34[i] + pow(T1(i,j),(-1-delta[i]))*pow((1-u(i,j)),theta[i])*log(1-u(i,j));
      //rowSums(T1(i,j)^delta[i])
      rowsum1[i] = rowsum1[i] + pow(T1(i,j),delta[i]);
      //rowSums(T1[, 2:1]*(T1+(1-u)^theta*(1+delta))*log(1-u))
      rowsum2[i] = rowsum2[i] + T1(i,2-j)*(T1(i,j) + pow((1-u(i,j)),theta[i])*(1+delta[i]))*log(1-u(i,j));
    }
  }

  for(i= 0;i<n;i++)
  {
    PT1(i,0) = T1(i,0)*T1(i,1);
    C1[i] = SD34[i]+pow(L1[i],(-(1+delta[i]/delta[i])))/(delta[i] - pow(L1[i],(-1/delta[i]))*delta[i]);
    C2[i] = (log(pow(L1[i],(1/delta[i]))-log(-1+pow(L1[i],(1/delta[i])))))/pow(theta[i],2);
  }

  for(i=0;i<n;i++)
  {
   outtheta[i] = (1/(pow(L1[i],3)*(-1-delta[i]*theta[i]+pow(L1[i],(1/delta[i]))*(1+delta[i])*theta[i]))*
      pow(PT1[i],(-1-2*delta[i]))*pow((rowsum1[i]-pow(PT1[i],delta[i])),2)*(
          1/theta[i]*pow(PT1[i],(-delta[i]))*(-SD12[i]*(1+delta[i])*theta[i]*(
              -2+theta[i]-2*delta[i]*theta[i]+pow(L1[i],(1/delta[i]))*(1+2*delta[i])*theta[i])-
              L1[i]*pow(PT1[i],(1+delta[i]))*(C1[i]*(-1+theta[i])*(
                  -1+theta[i]-theta[i]*delta[i]
                  +pow(L1[i],(1/delta[i]))*(1+delta[i])*theta[i])+
                  theta[i]*(C2[i]+delta[i]+C2[i]*delta[i]*theta[i]-
                  pow(L1[i],(1/delta[i]))*(1+delta[i])*(1+C2[i]*theta[i]))))+
                  L1[i]*(-1-delta[i]*theta[i]+pow(L1[i],(1/delta[i]))*(1+delta[i])*theta[i])
                     *rowsum2[i]));
  }


  return outtheta;

}
