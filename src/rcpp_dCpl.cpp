#include <Rcpp.h>
using namespace Rcpp;
NumericVector dCplrcpp(std::string CplNM, NumericMatrix u, List parCpl, bool log0);
NumericVector logDensFun(NumericMatrix u, NumericVector theta , NumericVector delta);
NumericVector dmvNormVecFun(int i,NumericVector x, NumericVector rho, Function f);
NumericVector dmvtRcpp(NumericVector x, NumericVector sigma, bool log0, Function f);
NumericVector vech2mRcpp(NumericVector vech, bool diag, Function f);

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
    int preBits;
    NumericVector density(u_nrow), out_logredoMPFR(u_nrow), out_logredo(u_nrow) ;
    delta = parCpl["delta"];
    theta = parCpl["theta"];

    // The usual log density
    out_log = logDensFun(u, theta, delta);
    // BB7 density is very unstable numerically. Use "Multiple Precision Floating-Point
    // Reliable" based on GNU Multiple Precision Library for "those errors only (NA, NAN,
    // Inf)" found in the result.
    Environment Rmpfr("package:Rmpfr");
    Function mpfr = Rmpfr["mpfr"];
    int redo_idx;
    int precBits = 1024;
    for(i=0;i<u_nrow;i++)
    {
      if(!std::isfinite(out_log[i]))
      {
        out_logredoMPFR[i] = logDensFun(u = mpfr(u(i,_), precBits = precBits),
                                        theta = mpfr(theta[i], precBits = precBits),
                                        delta = mpfr(delta[i], precBits = precBits));
        out_logredo[i] = out_logredoMPFR[i];
        //out_logredo[i] = Rcpp::as<float>(out_logredoMPFR[i]);
        out_log[i] = out_logredo[i];
      }
      else
      {
        warning("MPFR used with insufficient ", precBits, " precBits in BB7 density.");
      }
    }
  }

  if(CplNM == "gaussian")
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
    NumericVector logDensUpper(nObs),logDensLower(nObs), logDens(nObs) ;
    for(i=0;i<nObs;i++)
    {
      logDensUpper[i] = dmvNormVecFun(i,u_quantile(i,_),rho[i]);
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
    rho = parCpl["rho"]; //n-by-lq
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

    NumericVector logDensUpper(nObs),logDensLower(nObs), logDens(nObs) ;
    for(i=0;i<nObs;i++)
    {
      logDensUpper[i] = dmvt(u_quantile[i, , drop = FALSE],
                             vech2m(rho[i, ], diag = FALSE),
                             "shifted", // wikipedia type
                             df[i],TRUE);
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
