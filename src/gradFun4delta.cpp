#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gradFun4delta(NumericMatrix u, NumericVector theta , NumericVector delta)
{
  int i,j,n,m;
  n = u.nrow();
  m = u.ncol();

  NumericVector out(n),out1(n),out2(n), Tu1(n),Tv1(n),L1(n),L3(n),L4(n),D1(n),D2(n),L5(n);
  NumericMatrix T1(n,m), L34(n,m), D12(n,m);

  for(i= 0;i<n;i++)
  {
    for(j= 0;j<n;j++)
    {
      T1(i,j)  = 1-pow((1-u(i,j)),theta[i]);
      L34(i,j) = -1+pow(T1(i,j),delta[i]);
      D12(i,j) = pow(T1(i,j),(-2*delta[i]));
    }
  }

  for(i= 0;i<n;i++)
  {
    Tu1[i] = u(i,0);
    Tv1[i] = u(i,1);
    L1[i]  = -1 +pow(Tu1[i],delta[i]) + pow(Tv1[i],(-delta[i]));
    L3[i]  = L34(i,0);
    L4[i]  = L34(i,1);
    D1[i]  = D12(i,0);
    D2[i]  = D12(i,1);
    L5[i]  = pow(Tu1[i],delta[i])-L3[i]*pow(Tv1[i],delta[i]);
  }

  for(i=0;i<n;i++)
  {

    out1[i] =  (pow(L5[i],2)*D1[i]*D2[i]*
      ((-1+pow(L1[i],pow((1/delta[i]),2)))*pow(delta[i],2)*pow(theta[i],2)+log(L1[i])- 1/L5[i]*
      (-L5[i]*theta[i]*(delta[i]+pow(L1[i],(2/delta[i]))*(1+delta[i])*theta[i]-
      pow(L1[i],(1/delta[i]))*(3+delta[i]+(-1+delta[i])*theta[i]))*log(L1[i])

         + delta[i]*(pow(L1[i],(2/delta[i]))*(1+delta[i])*
           (L4[i]*pow(Tu1[i],delta[i])*delta[i]+pow(Tv1[i],delta[i])*(1+delta[i]))*pow(theta[i],2)
                       +(1+delta[i]*theta[i])*(pow(Tu1[i],delta[i])*delta[i]*theta[i]-pow(Tv1[i],delta[i])*(1+(1+pow(Tu1[i],delta[i]))*delta[i]*theta[i]))+
                         pow(L1[i],(1/delta[i]))*theta[i]*(L4[i]*pow(Tu1[i],delta[i])*delta[i]*(1+theta[i]+2*delta[i]*theta[i])
                                                             +pow(Tv1[i],delta[i])*(3-theta[i]+2*delta[i]*(1+theta[i]+theta[i]*delta[i]))))*log(Tu1[i])

         + delta[i]*(pow(L1[i],(2/delta[i]))*(1+delta[i])*
         (L3[i]*pow(Tv1[i],delta[i])*delta[i]+pow(Tu1[i],delta[i])*(1+delta[i]))*pow(theta[i],2)
                       +(1+delta[i]*theta[i])*(pow(Tv1[i],delta[i])*delta[i]*theta[i]-pow(Tv1[i],delta[i])*(1+(1+pow(Tv1[i],delta[i]))*delta[i]*theta[i]))+
                         pow(L1[i],(1/delta[i]))*theta[i]*(L3[i]*pow(Tv1[i],delta[i])*delta[i]*(1+theta[i]+2*delta[i]*theta[i])
                                                             +pow(Tu1[i],delta[i])*(3-theta[i]+2*delta[i]*(1+theta[i]+theta[i]*delta[i]))))*log(Tv1[i])

      )));



    out2[i] = ( pow(L1[i],2)*pow(delta[i],2)*theta[i]*(1+delta[i]*theta[i]+pow(L1[i],(2/delta[i]))*(1+delta[i])*theta[i]-
      pow(L1[i],(1/delta[i]))*(1+theta[i]+2*delta[i]*theta[i])));

    out[i]  =out1[i]/out2[i];

  }


  return out;

}
