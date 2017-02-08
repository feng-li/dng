#include <Rcpp.h>
using namespace Rcpp;
NumericVector dCpl(std::string CplNM, NumericMatrix u, List parCpl, bool log0);
NumericVector logDensFun(NumericMatrix u, NumericVector theta , NumericVector delta);
NumericVector dmvNormVecFun(int i,NumericVector x, NumericVector rho);
NumericVector dmvtRcpp(NumericVector x, NumericVector sigma, bool log0);

// [[Rcpp::export]]
NumericVector dCpl(std::string CplNM, NumericMatrix u, List parCpl, bool log0)
{
  Environment CalldCpl("package:dng");
  int i,j,n1,nObs,q;
  n1 = CplNM.size();
  nObs = u.nrow();
  q = u.ncol();
  NumericVector out_log(nObs), out(nObs);
  NumericVector delta(nObs),theta(nObs), df(nObs), u1(nObs) , u2(nObs);

  for(i=0;i<n1;i++)
  {if(CplNM[i]>=65 && CplNM[i]<=90)
    CplNM[i]=32+CplNM[i]; }

  if(CplNM == "bb7")
  {
    NumericVector density(nObs), out_logredoMPFR(nObs), out_logredo(nObs) ;
    delta = parCpl["delta"];
    theta = parCpl["theta"];

    // The usual log density
    out_log = logDensFun(u, theta, delta);
    // BB7 density is very unstable numerically. Use "Multiple Precision Floating-Point
    // Reliable" based on GNU Multiple Precision Library for "those errors only (NA, NAN,
    // Inf)" found in the result.
    Function logDensFun_infinite = CalldCpl["logDensFun_infinite"];
    NumericVector ldfinf(1);
    for(i=0;i<nObs;i++)
    {
      if(TRUE ==!R_FINITE(out_log[i]))
      {
        ldfinf  = logDensFun_infinite(u, theta, delta, i);
        out_log[i] = ldfinf[0];}
    }

  }

  else if(CplNM == "gaussian")
  {
    NumericMatrix rho = parCpl["rho"]; // n-by-lq
    NumericMatrix u_quantile(nObs,q); // n-by-q
    int lq = rho.ncol();
    NumericVector logDensUpper(nObs),logDensLower(nObs), logDens(nObs);
    NumericVector u_quantile_i(q),rhoi(lq);

    // The quantile for normal CDF
    for(i=0;i<nObs;i++)
    {
      for(j=0;j<q;j++)
      {
        u_quantile(i,j) = R::qnorm5(u(i,j), 0, 1, TRUE, FALSE);
        u_quantile_i[j] = u_quantile(i,j);
      }
    }
    for(i=0;i<nObs;i++)
    {
      for(j=0;j<lq;j++)
      { rhoi(j) = rho(i,j);}
    }

    //The CplNM density function C_12(u1, u2)
    NumericVector ldl;
    for(i=0;i<nObs;i++)
    {
      ldl = dmvNormVecFun(i,u_quantile_i,rhoi);
      logDensUpper[i] = ldl[0];
      logDensLower[i] = 0;
      for(j=0;j<q;j++)
      {logDensLower[i] = logDensLower[i] + R::dnorm4(u_quantile(i,j),0,1,TRUE);}
      logDens[i] = logDensUpper[i] - logDensLower[i];
    }
    //The output
    out_log = logDens;

  }

  else if(CplNM == "mtv")  //The multivariate t-copula
  { //Demarta & McNeil (2005),  The t copula and related copulas
    // df, corr
    df  = parCpl["df"];  //n-by-1
    NumericMatrix rho = parCpl["rho"]; //n-by-lq
    NumericMatrix u_quantile(nObs,q);
    int nObs;

    // The quantile for *univariate* t, that is x in t(x, df)
    for(i=0;i<nObs;i++)
    {
      for(j=0;j<q;j++)
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
    Environment callfun("package:dng");
    Function vech2m = callfun["vech2m"];
    NumericVector logDensUpper(nObs),logDensLower(nObs), logDens(nObs),u_quantile_i(q),rho_i(q),delta_vech2m(q);
    NumericVector uppv(1);
    for(i=0;i<nObs;i++)
    {
      for(j=0;j<q;i++)
      {
        u_quantile_i[j] = u_quantile(i,j);
        rho_i[j] = rho(i,j);
      }
      delta_vech2m = vech2m(rho_i, 0);
      uppv = dmvt(u_quantile_i,delta_vech2m, df[i],1, "shifted");
      logDensUpper[i] = uppv[0];

      logDensLower[i] = 0;
      for(j=0;j<q;j++)
      {logDensLower[i] = logDensLower[i] + R::dt(u_quantile(i,j),df[i],TRUE);}
      logDens[i] = logDensUpper[i] - logDensLower[i];
    }
    out_log = logDens;
  }
  else if(CplNM == "fgm")
  {
    NumericVector density(nObs);
    delta = parCpl["delta"];
    theta = parCpl["theta"];
    for(i=0;i<nObs;i++)
    {
      density[i] = 1+theta[i]*(1-2*u1[i])*(1-2*u2[i]);
      out[i] =  density[i];
    }

  }


  else if(CplNM == "gumbel")
  {
    NumericMatrix u_tilde(nObs,q);
    NumericVector u_tildeSumdelta(nObs), pctl_log(nObs);
    delta = parCpl["delta"];
    for(i=0;i<nObs ; i++)
    {
      for(j=0;j<q;j++)
      {
        u_tilde(i,j) = -log(u(i,j));
      }
    }

    for(i=0;i<nObs;i++)
    {
      u_tildeSumdelta[i] = 0;
      for(j=0;j<q;j++)
      {
        u_tildeSumdelta[i] = u_tildeSumdelta[i] + pow(u_tilde(i,j), delta[i]);
      }
    }
    for(i=0;i<nObs;i++)
    {
      pctl_log[i] = -pow(u_tildeSumdelta[i],1/delta[i]);
    }

    for(i=0; i<nObs;i++)
    {
      out_log[i] = (pctl_log[i]+u_tilde(i,0)+u_tilde(i,1)+
        (delta[i]-1)*(log(u_tilde(i,0))+log(u_tilde(i,1)))-
        (2-1/delta[i])*log(u_tildeSumdelta[i])+
        log(pow(u_tildeSumdelta[i],(1/delta[i]))+delta[i]-1));
    }


  }

  else if(CplNM == "clayton")
  {
    delta = parCpl["delta"];
    for(i=0; i<nObs; i++)
    {
      u1[i] = u(i,0);
      u2[i] = u(i,1);
    }
    for(i=0; i<nObs; i++)
    {
      out_log[i] = (log(1+delta[i]) + (-delta[i] -1)*(log(u1[i])+ log(u2[i])) +
        (-2-1/delta[i])*log( pow(u1[i],(-delta[i]))+pow(u2[i],(-delta[i]))-1));
    }

  }



  else
  {
    stop("Given copula name is not implemented.");
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
