#include <Rcpp.h>
using namespace Rcpp;
NumericVector dCplrcpp(std::string CplNM, NumericMatrix u, List parCpl, bool log0);
NumericVector logDensFun(NumericMatrix u, NumericVector theta , NumericVector delta);


// [[Rcpp::export]]
NumericVector dCpl_rcpp(std::string CplNM, NumericMatrix u, List parCpl, bool log0 = TRUE)
{
  int i,j,a,n1,n2,u_nrow,u_ncol;
  n1 = CplNM.size();
  u_nrow = u.nrow();
  u_ncol = u.ncol();
  NumericVector out_log(u_nrow), out(u_nrow);
  NumericVector delta(u_nrow),theta(u_nrow), df(u_nrow), rho(u_nrow), u1(u_nrow) , u2(u_nrow);

  for(i=0;i<n1;i++){
    if(CplNM[i]>=65 && CplNM[i]<=90)
      CplNM[i]=32+CplNM[i]; }


  if(CplNM == "bb7")
  {
    NumericVector density(u_nrow);
    delta = parCpl["delta"];
    theta = parCpl["theta"];

    // The usual log density
    out_log = logDensFun(u, theta, delta);

    // BB7 density is very unstable numerically. Use "Multiple Precision Floating-Point
    // Reliable" based on GNU Multiple Precision Library for "those errors only (NA, NAN,
    // Inf)" found in the result.


  }

  else if(CplNM == "gaussian")
  {
    int nObs;
    NumericMatrix u_quantile(u_nrow,u_ncol);
    rho = parCpl["rho"]; // n-by-lq
    // The quantile for normal CDF
    for(i=0;i<u_nrow;i++)
    {
      for(j=0;j<u_ncol;j++)
      {
        u_quantile(i,j) = R::qnorm5(u(i,j), 0, 1, TRUE, FALSE);
      }
    }
    nObs = u_quantile.nrow();

    //The CplNM density function C_12(u1, u2)
  }


  else if(CplNM == "mtv")
  {
    // df, corr
    df  = parCpl["df"];
    rho = parCpl["rho"];
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
    // u.quantile <- X

    // The log copula density function C_12(u1, u2)
    nObs = df.size();

    // the formula right before formula (1.4) in Genz and Bretz (2009), also in
    // Wikipedia.

    // The density of the t copula, Demarta & Department (2006) Eq(6)

    NumericMatrix logDensUpper(nObs,1),logDensLower(nObs,1), logDens(nObs,1) ;
    for(i=0;i<nObs;i++)
    {
      logDensUpper(i,0)=1;/* dmvt(x = u.quantile[i, , drop = FALSE],
                   sigma = vech2m(rho[i, ], diag = FALSE),
      type = "shifted", # wikipedia type
      df = df[i], log = TRUE)  */
    }

    //logDensLower <- apply(dt(u.quantile, df = df, log = TRUE), 1, sum)
    //logDens <- logDensUpper-logDensLower

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
