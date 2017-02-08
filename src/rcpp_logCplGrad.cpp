#include <Rcpp.h>
using namespace Rcpp;
List logCplGrad04(std::string CplNM, NumericMatrix u, List parCpl, std::string parCaller);
NumericVector gradFun4delta(NumericMatrix u, NumericVector theta , NumericVector delta);
NumericVector gradFun4theta(NumericMatrix u, NumericVector theta , NumericVector delta);
NumericVector gradFun4df(int i, NumericMatrix rho, NumericVector df , NumericMatrix u_quantile);
List MargiModelGrad(NumericVector y, List par,std::string type, std::string parCaller,GenericVector denscaller);
NumericMatrix FUN(int i, NumericMatrix x, NumericVector mu, NumericVector df, NumericMatrix rho, NumericVector uIdx );
// [[Rcpp::export]]
List logCplGrad(std::string CplNM, NumericMatrix u, List parCpl, std::string parCaller)
{
  Environment Callfunc("package:logCplGradcpp");
  int i,j;
  int n_parCaller = parCaller.size();
  int n_cplnm = CplNM.size();
  int nObs = u.nrow();
  int q    = u.ncol();
  //char parcaller[n_parCaller];
  for(i=0;i<n_parCaller;i++){
    if(parCaller[i]>=65 && parCaller[i]<=90)
      parCaller[i]=32+parCaller[i]; }
  const char *parcal = parCaller.c_str();
  for(i=0;i<n_cplnm;i++){
    if(CplNM[i]>=65 && CplNM[i]<=90)
      CplNM[i]=32+CplNM[i]; }

  NumericVector out_delta, out_theta, out_df, out_u;
  NumericMatrix out_rho;


  if(CplNM == "bb7")
  {
    // The standard copula parameters (recycled if necessary, should be a vector).
    NumericVector    delta = parCpl["delta"];
    NumericVector    theta = parCpl["theta"];
    //NumericVector A,gradout_redoMPFR,gradout_redo;

    if(strstr(parcal,"delta"))
    {
      //################################################################################
      //### DEBUGGING
      //## u <- matrix(c(0.6, 0.3), 1, )
      //## theta <- 3.5
      //## delta <- 2.4
      //### PASSED
      //################################################################################
      NumericVector dinf(1);
      NumericVector gradout = gradFun4delta(u = u, theta = theta, delta = delta);
      Function gradFun4delta_infinite = Callfunc["gradFun4delta_infinite"];
      for(i=0;i<nObs;i++)
      {
        if(TRUE ==!R_FINITE(gradout[i]))
        {
          dinf=gradFun4delta_infinite(u, delta, theta, i, parCaller);
          gradout[i]  = dinf[0];
        }

      }
      out_delta = gradout;

    }

    if(strstr(parcal,"theta"))
    {
      //################################################################################
      //### DEBUGGING
      //## u <- matrix(c(0.6, 0.3), 1, )
      //## theta <- 3.5
      //## delta <- 2.4
      //### PASSED
      //################################################################################
      NumericVector gradout = gradFun4theta(u = u, theta = theta, delta = delta);
      Function gradFun4theta_infinite = Callfunc["gradFun4theta_infinite"];
      NumericVector tinf(1);
      for(i=0;i<nObs;i++)
      {
        if(TRUE ==!R_FINITE(gradout[i]))
        {
          tinf  = gradFun4theta_infinite(u, delta, theta, i, parCaller);
          gradout[i] = tinf[0];
        }
      }
      out_theta = gradout;
    }

    if(strstr(parcal,"u"))
    {
      //## Gradient w.r.t u. NOTE: The BB7 copula's marginal are exchangeable which means
      //## the expression for the gradient w.r.t u1 and u2 are the same if swap u1 and u2

      //################################################################################
      //## DEBUGGING
      //## u <- matrix(c(0.2, 0.3), 1, )
      //## theta <- 1.5
      //## delta <- 2.4
      //## PASSED
      //################################################################################
      NumericMatrix ub(nObs,q),D12(nObs,q);
      NumericVector gradCpl_u(nObs), ub1(nObs), D1(nObs), D2(nObs), S1(nObs), S2(nObs), S3(nObs);
      for(i=0;i<nObs;i++)
      {
        S1[i] = 1;
        S3[i] = 0;
        for(j=0;j<q;j++)
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
      out_u = gradCpl_u;

    }

  }

  else if(CplNM == "mvt")
  {
    NumericVector df = parCpl["df"];
    NumericMatrix rho = parCpl["rho"];

    NumericMatrix u_quantile(nObs,q);
    NumericVector u_quantile_j(q);
    for(j=0; j<q; j++)
    {
      for(i=0; i<nObs; i++)
      { u_quantile(i,j) = R::qt(u(i,j),df(i), TRUE, FALSE);
        u_quantile_j(i) = u_quantile(i,j);}
    }


    if(strstr(parcal,"df"))
    {
      NumericVector logCplGrad_df(nObs), logCplGrad_df_upper(nObs),logCplGrad_df_lower(nObs),logCplGrad_df_lowerMat(q);
      NumericVector dfupp(1);
      List mat;
      List parM = List::create(_["mu"] = 0,
                               _["df"] = df,
                               _["phi"] = 1,
                               _["lmd"] = 1);
      for(i=0; i<nObs; i++)
      {
        dfupp = gradFun4df(i,rho,df,u_quantile) ;
        logCplGrad_df_upper[i] = dfupp[0];
        logCplGrad_df_lower[i] = 0;
        for(j=0; j<q;j++)
        {
          mat = MargiModelGrad( u_quantile_j,parM , "splitt", "df", "d");
          logCplGrad_df_lowerMat = mat["d"];
          logCplGrad_df_lower[i] = logCplGrad_df_lower[i]+logCplGrad_df_lowerMat[j];
        }
        logCplGrad_df[i] = logCplGrad_df_upper[i] -logCplGrad_df_lower[i];
      }
      out_df = logCplGrad_df;
    }

    if(strstr(parcal,"rho"))
    {
      Environment Callfunc("package:grand");
      Function mult = Callfunc["mult"];
      Function vech2m = Callfunc["vech2m"];
      Environment base("package:base");
      Function solve = base["solve"];
      Function t = base["t"];

      int lq = rho.ncol();
      NumericMatrix logCplGrad_rho(nObs, lq);  // n-by-lq
      NumericMatrix Sigma;
      NumericVector rhoi(lq);
      for(j=0;j<lq;j++)
      {    rhoi[j] = rho(i,j);  }
      for(i=0;i<nObs;i++)
      {
        Sigma = vech2m(rhoi, 0);
        double v = df[i];
        NumericMatrix x(1,q);  //col-vector
        for(j=0;j<q;j++)
        { x(0,j) = u_quantile(i,j);  }
        int mu = 0;
        int p = Sigma.nrow();
        int q = Sigma.ncol();

        // C0 <- as.vector(t(x-mu)%*%solve(Sigma)%*%(x-mu))
        NumericMatrix C1 = solve(Sigma, (x-mu));
        NumericVector C0 = mult(t(x-mu),C1);
        NumericMatrix solveSigma = solve(Sigma);
        NumericMatrix multC1tC1 = mult((-C1),t(C1));
        NumericMatrix logGradCpl_Sigma(p,q);
        for(int a;a<p;a++)
        {
          for(int b;b<q;b++)
          {
            logGradCpl_Sigma(a,b) = (-1/2*solveSigma(a,b)-(v+p)/2*pow((1+C0[b]/v),(-1))*multC1tC1(a,b)/v);
          }
        }

        for(int a=0;a<lq;a++)
        {
          for(int b= a+1;b<lq;b++)
          {logCplGrad_rho(a,b) = logGradCpl_Sigma(a,b);}
        }
    out_rho = logCplGrad_rho;

    }


    }

    if(strstr(parcal,"u"))
    {
      NumericVector logCplGrad_u(nObs);
      // The gradient with respect to u_i Reorder the parameters.
      //double imar = substr(parCaller, 2, nchar(parCaller));
      double imar = 0;
      NumericVector uIdx(q) ;
      for(i=0;i<q;i++)
      { uIdx[i] = i+1; }
      uIdx[0] = imar;
      uIdx[imar] = 1;

      //u <- u[, uIdx]
      //x <- u.quantile[, uIdx]
      //x1 <- x[, 1]
      NumericMatrix x;
      NumericVector x1;
      for(int i=0;i<nObs;i++)
      {
        for(j=0;j<q;j++)
        {
          x(i,j) = u_quantile(i,uIdx[j]);
        }
        x1[i] = x(i,1);
      }

      int mu = 0;
      // The t copula and related copulas EQ.(6)
      NumericMatrix fx1(nObs,1) ; // n-by-1
      NumericMatrix f1x1(nObs,q); // n-by-q
      for(i=0;i<nObs;i++)
      {
        fx1(i,0) = R::dt(x1[i],df[i],FALSE);
        for(int j=0;j<q;j++)
        {
          f1x1(i,j) = -(x1[i]-mu)*df[i]*(1+df[i])*pow((df[i]/(df[i]+pow((x1[i]-mu),2))),(-1+(1+df[i]/2)))/
            (pow((df[i]+pow((x1[i]-mu),2)),2)*sqrt(df[i])*R::beta(df[i]/2,1/2));
        }
      }

      // The first marginal CDF derivative with respect to x1.
      NumericMatrix F1x1(nObs,1) ; //n-by-1
      NumericVector I0(nObs),I(nObs);
      for(i=0;i<nObs;i++)
      {
        I0[i]=(x1[i]<mu);
        I[i]= 1-I0[i];
      }
      for(i=0;i<nObs;i++)
      {
        if(x1[i]>=mu)
        {
          F1x1(i,0) =  (-(2*pow((x1[i]-mu),3)/(pow((x1[i]-mu),2)+pow(df[i],2)-2*(x1[i]-mu)/(pow((x1[i]-mu),2)+df[i]))*
            pow((1-pow((x1[i]-mu),2)/(pow((x1[i]-mu),2)+df[i])),(-1+df[i]/2))/
              (2*sqrt(pow((x1[i]-mu),2)/(pow((x1[i]-mu),2)+df[i]))*
                R::beta(1/2, df[i]/2))));
        }
        else
        {
          F1x1(i,0) = ( -(x1[i]-mu)*df[i]*
            pow((df[i]/(pow((x1[i]-mu),2)+df[i])),(-1+df[i]/2))/pow((pow((x1[i]-mu),2)+df[i]),2)/
              sqrt(1-df[i]/(pow((x1[i]-mu),2)+df[i]))/R::beta(df[i]/2, 1/2));
        }
      }

      NumericMatrix gradLogCpl_x1(nObs,1);
      NumericMatrix funM;
      for(i=0;i<nObs;i++)
      {
        funM = FUN(i, x, mu, df, rho, uIdx);
        gradLogCpl_x1(i,0) = funM(0,0);
      }
      logCplGrad_u = gradLogCpl_x1*(1/F1x1)- 1/fx1*f1x1/F1x1;
      out_u = logCplGrad_u;
    }

  }

  else if(CplNM == "gumbel")
  {
    NumericVector A(nObs),u1t(nObs),u2t(nObs),uDelta(nObs);
    NumericVector A1(nObs), A2(nObs), A3(nObs), A4(nObs), A5(nObs), A6(nObs);

    NumericVector delta = parCpl["delta"];
    for(i=0; i<nObs; i++)
    {
      u1t[i] = -log(u(i,0));
      u2t[i] = -log(u(i,1));
      uDelta[i] = pow(u1t[i],delta[i]) + pow(u2t[i],delta[i]);
    }

    if(strstr(parcal,"delta"))
    {
      NumericVector logCplGrad_delta(nObs),u1(nObs),u2(nObs);
      for(i=0; i<nObs; i++)
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
      out_delta = logCplGrad_delta;
    }

    if(strstr(parcal,"u"))
    {
      NumericVector logCplGrad_u(nObs), u1(nObs);
      for(i=0; i<nObs; i++)
      {
        u1[i] = u(i,0);
        logCplGrad_u[i] = (-1/u1[i] + ((-1+delta[i] - (1+pow(uDelta[i],(2/delta[i]))
                                                         + 3*pow(uDelta[i],(1/delta[i]))*(-1+delta[i])
                                                         - 3*delta[i]+2*pow(delta[i],2))*pow((u1t[i]),delta[i])/(
                                                             uDelta[i]*(-1+pow(uDelta[i],(1/delta[i]))+delta[i])))/(-u1t[i]*u1[i])));

      }
      out_u = logCplGrad_u;
    }

    }
  else if(CplNM == "clayton")
  {
    NumericVector delta = parCpl["delta"];
    NumericVector u1(nObs), u2(nObs);
    for(i=0; i<nObs; i++)
    {
      u1[i] = u(i,0);
      u2[i] = u(i,1);
    }

    if(strstr(parcal,"delta"))
    {
      NumericVector A(nObs),logCplGrad_delta(nObs);
      for(i=0; i<nObs; i++)
      {
        A[i] = (-1 + pow(u1[i],(-delta[i])) + pow(u2[i],(-delta[i])));

        logCplGrad_delta[i] = (1/(1 + delta[i]) - log(u1[i]) - log(u2[i]) +
          (-2 - 1/delta[i])*(-pow(u1[i],(-delta[i]))*log(u1[i])
                               - pow(u2[i],(-delta[i]))*log(u2[i]))/A[i] +
                                 log(A[i])/pow(delta[i],2));
      }
      out_delta = logCplGrad_delta;
    }
    if(strstr(parcal,"u"))
    {
      NumericVector A(nObs),logCplGrad_u(nObs);
      for(i=0; i<nObs;i++)
      {
        A[i] = pow(u1[i],delta[i])*(-1+pow(u2[i],delta[i]));

        logCplGrad_u[i] = (-(pow(u2[i],delta[i])*delta[i] + A[i]*(1+delta[i]))/(u1[i]*(-pow(u2[i],delta[i])+ A[i])));
      }
      out_u = logCplGrad_u;
    }
    }
  else
  {
    stop("No such copula defined!");
  }


  List out = List::create(_["delta"] = out_delta,
                          _["theta"] = out_theta,
                          _["df"] = out_df,
                          _["rho"] = out_rho,
                          _["u"] = out_u);

  return out;


}
