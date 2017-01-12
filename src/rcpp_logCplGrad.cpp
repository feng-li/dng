#include <Rcpp.h>
using namespace Rcpp;
List logCplGrad(std::string CplNM, NumericMatrix u, List parCpl, std::string parCaller);
NumericVector gradFun4delta(NumericMatrix u, NumericVector theta , NumericVector delta);
NumericVector gradFun4theta(NumericMatrix u, NumericVector theta , NumericVector delta);
NumericVector gradFun4df(int i, NumericMatrix rho, NumericVector df , NumericMatrix u_quantile);

// [[Rcpp::export]]

List logCplGrad(std::string CplNM, NumericMatrix u, List parCpl, std::string parCaller)
{
  Environment Callfunc("package:logCplGrad02");
  int i,j,a;
  int n_parCaller = parCaller.size();
  int n_cplnm = CplNM.size();
  int u_nrow = u.nrow();
  int u_ncol = u.ncol();
  int q = u.ncol();
  NumericVector out_log(u_nrow), logCplGrad_delta(u_nrow), logCplGrad_u(u_nrow),logCplGrad_df(u_nrow)£¬logCplGrad_rho(u_nrow);
  NumericVector delta(u_nrow),theta(u_nrow), gradout(u_nrow), df(u_nrow),rho(u_nrow), u1(u_nrow) , u2(u_nrow),gradCpl_u(u_nrow);


  for(i=0;i<n_parCaller;i++){
    if(parCaller[i]>=65 && parCaller[i]<=90)
      parCaller[i]=32+parCaller[i]; }

  for(i=0;i<n_cplnm;i++){
    if(CplNM[i]>=65 && CplNM[i]<=90)
      CplNM[i]=32+CplNM[i]; }


  if(CplNM == "bb7")
  {
    NumericVector A,gradout_redoMPFR,gradout_redo;
    int precBits;
    delta = parCpl["delta"];
    theta = parCpl["theta"];

//if( "delta" %in% tolower(parCaller))
    //################################################################################
    //### DEBUGGING
    //## u <- matrix(c(0.6, 0.3), 1, )
    //## theta <- 3.5
    //## delta <- 2.4
    //### PASSED
    //################################################################################
    gradout = gradFun4delta(u = u, theta = theta, delta = delta);
    Function gradFun4delta_infinite = Callfunc["gradFun4delta_infinite"];
    for(i=0;i<u_nrow;i++)
    {
      if(TRUE ==!R_FINITE(gradout[i]))
      { gradout[i]  = gradFun4delta_infinite(u, delta, theta, i);}
    }

//if( "theta" %in% tolower(parCaller))
      //################################################################################
      //### DEBUGGING
      //## u <- matrix(c(0.6, 0.3), 1, )
      //## theta <- 3.5
      //## delta <- 2.4
      //### PASSED
      //################################################################################
    gradout = gradFun4theta(u = u, theta = theta, delta = delta);
    Function gradFun4theta_infinite = Callfunc["gradFun4theta_infinite"];
    for(i=0;i<u_nrow;i++)
    {
      if(TRUE ==!R_FINITE(gradout[i]))
      { gradout[i]  = gradFun4theta_infinite(u, delta, theta, i);}
    }
    



    //## Gradient w.r.t u. NOTE: The BB7 copula's marginal are exchangeable which means
    //## the expression for the gradient w.r.t u1 and u2 are the same if swap u1 and u2

    //################################################################################
    //## DEBUGGING
    //## u <- matrix(c(0.2, 0.3), 1, )
    //## theta <- 1.5
    //## delta <- 2.4
    //## PASSED
    //################################################################################
    NumericMatrix ub(u_nrow,u_ncol),D12(u_nrow,u_ncol);
    NumericVector gradCpl_u(u_nrow), ub1(u_nrow), D1(u_nrow), D2(u_nrow), S1(u_nrow), S2(u_nrow), S3(u_nrow);
    for(i=0;i<u_nrow;i++)
    {
      S1[i] = 1;
      S3[i] = 0;
      for(j=0;j<u_ncol;j++)
      {
        ub(i,j)  = 1 - u(i,j);
        D12(i,j) = 1 - pow(ub(i,j),theta[i]);
        S1[i] = S1[i] - pow(D12(i,j),(-delta[i]));
        S2[i] = 1 - pow((-S1[i]),(1/delta[i]));
        S3[i] = S3[i] + pow(D12(i,j),delta[i]);
      }

      ub1[i] = ub(i,0);
      D1[i] =  D12(i,0);
      D2[i] =  D12(i,1);
      gradCpl_u[i] = (
        -(pow(D1[i],(-1-3*delta[i]))*pow(D2[i],(-2*delta[i]))*
          pow((S3[i]-pow(D1[i],delta[i])*pow(D2[i],delta[i])),2)*
          (pow(D1[i],(1+delta[i]))*S1[i]*(-1+theta[i])*
          (1+pow((-S1[i]),(1/delta[i]))*(-1+theta[i])-
          theta[i]+pow(S2[i],2)*theta[i]+pow(S2[i],2)*delta[i]*theta[i])+
          pow(D2[i],(-delta[i]))*(pow(D1[i],delta[i])*(-1+pow(D2[i],delta[i]))*S2[i]*(1+delta[i])*theta[i]*
          (-1+(1+S2[i]+S2[i]*delta[i])*theta[i])+pow(D2[i],delta[i])*
          (1+theta[i]*(-3+2*theta[i]+S2[i]*(1+delta[i])*
          (-2+(2+S2[i]*delta[i])*theta[i]))))*pow(ub1[i],theta[i])))/
      (pow(S1[i],3)*(1+delta[i]*theta[i]+pow((-S1[i]),(2/delta[i]))*(1+delta[i])*theta[i]-
        pow((-S1[i]),(1/delta[i]))*(1+theta[i]+2*theta[i]*delta[i]))*ub1[i])
      );
    }

}

  else if(CplNM == "mvt")
  {
    int nObs;
    nObs = u.nrow();
    NumericVector logCplGrad_df_upper(u_nrow), logCplGrad_df_lowerMat(u_nrow);
    NumericMatrix u_quantile(u_nrow,u_ncol);
    df = parCpl["rho"];
    rho = parCpl["rho"];
    for(j=0; j<u_ncol; j++)
    {
      for(i=0; i<u_nrow; i++)
      { u_quantile[i] = R::qt(u(i,j),df(i), TRUE, FALSE);}
    }

    for(i=0; i<u_nrow; i++)
    {
      logCplGrad_df_upper[i] = 0 ;
      logCplGrad_df_lowerMat[i] = 0 ;
    }

  }

  else if(CplNM == "gumbel")
  {
    NumericVector A(u_nrow),u1t(u_nrow),u2t(u_nrow),uDelta(u_nrow);
    NumericVector A1(u_nrow), A2(u_nrow), A3(u_nrow), A4(u_nrow), A5(u_nrow), A6(u_nrow);

    delta = parCpl["delta"];
    for(i=0; i<u_nrow; i++)
    {
      u1t[i] = -log(u(i,0));
      u2t[i] = -log(u(i,1));
      uDelta[i] = pow(u1t[i],delta[i]) + pow(u2t[i],delta[i]);
    }

    for(i=0; i<u_nrow; i++)
    {
      u1[i] = u(i,0);
      u2[i] = u(i,1);

      A1[i] = (-1+pow(uDelta[i],(1/delta[i]))+delta[i]);
      A2[i] = uDelta[i]*log(uDelta[i]);
      A3[i] = pow(u1t[i],delta[i])*log(u1t[i]);
      A4[i] = pow(u2t[i],delta[i])*log(u2t[i]);
      A6[i] = pow(uDelta[i],(-1+1/delta[i]))*(A2[i] - (A3[i] + A4[i])*delta[i])/pow(delta[i],2);


      logCplGrad_delta[i] = (-log(uDelta[i])/pow(delta[i],2)+log(u1t[i]) + log(u2t[i])
                               + (1-2*delta[i])*(A3[i]+A4[i])/(uDelta[i]*delta[i]) + A6[i] + (1 -A6[i])/A1[i]);
    }

    for(i=0; i<u_nrow; i++)
    {
      u1[i] = u(i,0);
      logCplGrad_u[i] = (-1/u1[i] + ((-1+delta[i] - (1+pow(uDelta[i],(2/delta[i]))
                                                       + 3*pow(uDelta[i],(1/delta[i]))*(-1+delta[i])
                                                       - 3*delta[i]+2*pow(delta[i],2))*pow((u1t[i]),delta[i])/(
                                                           uDelta[i]*(-1+pow(uDelta[i],(1/delta[i]))+delta[i])))/(-u1t[i]*u1[i])));

    }

  }


  else if(CplNM == "clayton")
  {
    NumericVector A;
    delta = parCpl["delta"];
    for(i=0; i<u_nrow; i++)
    {
      u1[i] = u(i,0);
      u2[i] = u(i,1);
    }
    for(i=0; i<u_nrow; i++)
    {
      A[i] = (-1 + pow(u1[i],(-delta[i])) + pow(u2[i],(-delta[i])));

      logCplGrad_delta[i] = (1/(1 + delta[i]) - log(u1[i]) - log(u2[i]) +
        (-2 - 1/delta[i])*(-pow(u1[i],(-delta[i]))*log(u1[i])
                             - pow(u2[i],(-delta[i]))*log(u2[i]))/A[i] +
                               log(A[i])/pow(delta[i],2));
    }

    for(i=0; i<u_nrow;i++)
    {
      A[i] = pow(u1[i],delta[i])*(-1+pow(u2[i],delta[i]));

      logCplGrad_u[i] = (-(pow(u2[i],delta[i])*delta[i] + A[i]*(1+delta[i]))/(u1[i]*(-pow(u2[i],delta[i])+ A[i])));

    }

  }

  List out = List::create(_["delta"] = gradout,
                          _["theta"] = gradout,
                          _["u"] = gradCpl_u,
                          _["df"] = logCplGrad_df,
                          _["rho"] =  logCplGrad_rho
                            );

  return out;
}
