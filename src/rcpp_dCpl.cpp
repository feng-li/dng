#include <Rcpp.h>
using namespace Rcpp;
NumericVector dCplrcpp(std::string CplNM, NumericMatrix u, List parCpl, bool log0);
NumericVector logDensFun(NumericMatrix u, NumericVector theta , NumericVector delta);
double dmvNormVecFun(int i,NumericVector x, NumericVector rho);
NumericVector dmvtRcpp(NumericVector x, NumericVector sigma, bool log0);
NumericVector vech2mCpp(NumericVector vech, bool diag);

// [[Rcpp::export]]
NumericVector dCplrcpp(std::string CplNM, NumericMatrix u, List parCpl, bool log0)
{
  Environment CalldCpl("package:dCpl002");
  int i,j,n1,u_nrow,u_ncol;
  n1 = CplNM.size();
  u_nrow = u.nrow();
  u_ncol = u.ncol();
  NumericVector out_log(u_nrow), out(u_nrow);
  NumericVector delta(u_nrow),theta(u_nrow), df(u_nrow), u1(u_nrow) , u2(u_nrow);


  for(i=0;i<n1;i++){
    if(CplNM[i]>=65 && CplNM[i]<=90)
      CplNM[i]=32+CplNM[i]; }


  if(CplNM == "bb7")
  {

    NumericVector density(u_nrow), out_logredoMPFR(u_nrow), out_logredo(u_nrow) ;
    delta = parCpl["delta"];
    theta = parCpl["theta"];

    // The usual log density
    out_log = logDensFun(u, theta, delta);
    // BB7 density is very unstable numerically. Use "Multiple Precision Floating-Point
    // Reliable" based on GNU Multiple Precision Library for "those errors only (NA, NAN,
    // Inf)" found in the result.
    Function logDensFun_infinite = CalldCpl["logDensFun_infinite"];
    NumericVector ldfinf(1);
    for(i=0;i<u_nrow;i++)
    {
      if(TRUE ==!R_FINITE(out_log[i]))
      {
       ldfinf  = logDensFun_infinite(u, theta, delta, i);
       out_log[i] = ldfinf[0];}
    }

  }

  else if(CplNM == "gaussian")
  {
    int nObs,lq;
    NumericMatrix rho = parCpl["rho"]; // n-by-lq
    NumericMatrix u_quantile(u_nrow,u_ncol); // n-by-q
    NumericVector logDensUpper(u_nrow),logDensLower(u_nrow), logDens(u_nrow);
    NumericVector u_quantile_i(u_ncol),rhoi(lq);
    // The quantile for normal CDF
    for(i=0;i<u_nrow;i++)
    {
      for(j=0;j<u_ncol;j++)
      {
        u_quantile(i,j) = R::qnorm5(u(i,j), 0, 1, TRUE, FALSE);
        u_quantile_i = u_quantile(i,j);
      }
    }
    nObs = u_quantile.nrow();
    lq = rho.ncol();
    for(i=0;i<u_nrow;i++)
    {
      for(j=0;j<lq;j++)
      {
        rhoi(j) = rho(i,j);
      }
    }
    //The CplNM density function C_12(u1, u2)
    for(i=0;i<nObs;i++)
    {
      logDensUpper[i] = dmvNormVecFun(i,u_quantile_i,rhoi);
      logDensLower[i] = 0;
      for(j=0;j<u_ncol;j++)
      {logDensLower[i] = logDensLower[i] + R::dnorm4(u_quantile(i,j),0,1,TRUE);}
      logDens[i] = logDensUpper[i] - logDensLower[i];
      //The output
      out_log = Rcpp::as<NumericMatrix>(logDens);
    }
  }

  else if(CplNM == "mtv")  //The multivariate t-copula
  { //Demarta & McNeil (2005),  The t copula and related copulas
    // df, corr
    df  = parCpl["df"];  //n-by-1
    NumericMatrix rho = parCpl["rho"]; //n-by-lq
    NumericMatrix u_quantile(u_nrow,u_ncol);
    int nObs;

    // The quantile for *univariate* t, that is x in t(x, df)
    for(i=0;i<u_nrow;i++)
    {
      for(j=0;j<u_ncol;j++)
      {
        u_quantile(i,j) = R::qt(u(i,j),df[i],TRUE,FALSE);
      }
    }
    // u_quantile = X

    // The log copula density function C_12(u1, u2)
    nObs = df.size();

    // the formula right before formula (1.4) in Genz and Bretz (2009), also in
    // Wikipedia.

    // The density of the t copula, Demarta & Department (2006) Eq(6)

    Environment mvtnorm("package:mvtnorm");
    Function dmvt = mvtnorm["dmvt"];
    NumericVector logDensUpper(nObs),logDensLower(nObs), logDens(nObs),u_quantile_i(u_ncol),rho_i(u_ncol),delta_vech2m(u_ncol);
    NumericVector uppv(1);
    for(i=0;i<nObs;i++)
    {
      for(j=0;j<u_ncol;i++)
      {
        u_quantile_i[j] = u_quantile(i,j);
        rho_i[j] = rho(i,j);
        }
      delta_vech2m = vech2mCpp(rho_i, FALSE);
      uppv = dmvt(u_quantile_i,delta_vech2m, df[i],1, "shifted");
      logDensUpper[i] = uppv[0];

      logDensLower[i] = 0;
      for(j=0;j<u_ncol;j++)
      {logDensLower[i] = logDensLower[i] + R::dt(u_quantile(i,j),df[i],TRUE);}
      logDens[i] = logDensUpper[i] - logDensLower[i];
    }

  }

  else if(CplNM == "fgm")
  {
    NumericVector density(u_nrow);
    delta = parCpl["delta"];
    theta = parCpl["theta"];
    for(i=0;i<u_nrow;i++)
    {
      density[i] = 1+theta[i]*(1-2*u1[i])*(1-2*u2[i]);
      out[i] =  density[i];
    }

  }

  else if(CplNM == "gumbel")
  {
    NumericMatrix u_tilde(u_nrow,u_ncol);
    NumericVector u_tildeSumdelta(u_nrow), pctl_log(u_nrow);
    delta = parCpl["delta"];
    for(i=0;i<u_nrow ; i++)
    {
      for(j=0;j<u_ncol;j++)
      {
        u_tilde(i,j) = -log(u(i,j));
      }
    }

    for(i=0;i<u_nrow;i++)
    {
      u_tildeSumdelta[i] = 0;
      for(j=0;j<u_ncol;j++)
      {
        u_tildeSumdelta[i] = u_tildeSumdelta[i] + pow(u_tilde(i,j), delta[i]);
      }
    }
    for(i=0;i<u_nrow;i++)
    {
      pctl_log[i] = -pow(u_tildeSumdelta[i],1/delta[i]);
    }

    for(i=0; i<u_nrow;i++)
    {
      out_log[i] = (pctl_log[i]+u_tilde(i,0)+u_tilde(i,1)+
        (delta[i]-1)*(log(u_tilde(i,0))+log(u_tilde(i,1)))-
        (2-1/delta[i])*log(u_tildeSumdelta[i])+
        log(pow(u_tildeSumdelta[i],(1/delta[i]))+delta[i]-1));
    }


  }


  else if(CplNM == "frank")
  {
    NumericVector eta(u_ncol);
    delta = parCpl["delta"];
    for(i=0; i<u_nrow; i++)
    {
      u1[i] = u(i,0);
      u2[i] = u(i,1);
    }
    for(i=0; i<u_nrow; i++)
    {
      eta[i] = 1-exp(-delta[i]);
      out_log[i] = (log(delta[i])+log(eta[i])-delta[i]*(u1[i]+u2[i])-
        2*log(eta[i]-(1-exp(-delta[i]*u1[i]))*-(1-exp(-delta[i]*u1[i]))));
    }
  }


  else if(CplNM == "clayton")
  {
    delta = parCpl["delta"];
    for(i=0; i<u_nrow; i++)
    {
      u1[i] = u(i,0);
      u2[i] = u(i,1);
    }
    for(i=0; i<u_nrow; i++)
    {
      out_log[i] = (log(1+delta[i]) + (-delta[i] -1)*(log(u1[i])+ log(u2[i])) +
        (-2-1/delta[i])*log( pow(u1[i],(-delta[i]))+pow(u2[i],(-delta[i]))-1));
    }

  }


  if(log0)
  {
    out = out_log;
  }
  else
  {
    out = exp(out_log);
  }

  return out;
}
