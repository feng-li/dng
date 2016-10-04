#include <Rcpp.h>
using namespace Rcpp;
List MargiModelGrad_rcpp(NumericVector y, List par,std::string type, std::string parCaller,GenericVector denscaller);
NumericMatrix ghypergeo(NumericMatrix a, NumericMatrix b, NumericVector z,int k);
double ibeta(double x, double a, double b, bool log0, bool reg);

// [[Rcpp::export]]
List MargiModelGrad_rcpp(NumericVector y, List par,std::string type, std::string parCaller,GenericVector denscaller)
{
  int i,j,a,n1,n2,n3,n4;

  n1=y.size();
  n2=type.size();
  n3=parCaller.size();
  n4=denscaller.size();

  NumericVector mu(n1),df(n1),phi(n1),lmd(n1),logMargiDens(n1),outu(n1),outd(n1);


  for(i=0;i<n2;i++){
    if(parCaller[i]>=65 && type[i]<=90)
      parCaller[i]=32+type[i]; }

  for(i=0;i<n3;i++){
    if(parCaller[i]>=65 && parCaller[i]<=90)
      parCaller[i]=32+parCaller[i]; }



  ////////////////////////////gaussian////////////////////////////////////////
  if(type == "gaussian")
  {
    mu  = par["mu"];
    phi = par["phi"];
    for(i=0;i<n1;i++)
    {logMargiDens[i]=R::dnorm4(y[i],mu[i],phi[i], TRUE);}

    if(parCaller == "mu")
    {
      for(j=0;j<n1;j++)
      {
        outu[j] = exp(logMargiDens[j]);
        outd[j] = (y[j-1]-mu[j])/pow(phi[j],2);
      }

    }

    else if(parCaller == "phi")
    {
      for(j=0;j<n1;j++)
      {
        outu[j] = -(y[j]-mu[j])/phi[j]*exp(logMargiDens[j]);
        outd[j] = pow((y[j]-mu[j]),2)/pow(phi[j],3)-1/phi[j];
      }

    }

    else
    { stop("No such parameter!");}

  }

  /////////////////////////////splitt///////////////////////////////////////

  else if(type == "splitt")
  {
    mu  = par["mu"];
    df  = par["df"];
    phi = par["phi"];
    lmd = par["lmd"];
    // PDF
    // NumericVector MargiDens;
    // MargiDens[i] = R::dsplitt(y, mu,  df,  phi,  lmd, FALSE)

    if(parCaller == "mu")
    {
      if(std::strncmp(denscaller[0],"u",1)==0)
      {
        NumericVector I0(n1),I(n1), sign(n1),beta0(n1);

        for(int a=0;a<n1;a++)
        {
          I0[a]=(y[a]<=mu[a]);
          I[a]= 1-I0[a];
          sign[a]=1*I0[a]+lmd[a]*I[a];
          beta0[a]=R::beta(df[a]/2, 1/2);
        }

        for(j=0;j<n1;j++)
        {
          outu[j] = -2*sign[j]*sqrt(1/(pow((y[j]-mu[j]),2)+pow(sign[j],2)*df[j]*pow(phi[j],2)))*
            pow((pow(sign[j],2)*df[j]*pow(phi[j],2)/
                   (pow((y[j]-mu[j]),2)+pow(sign[j],2)*df[j]*pow(phi[j],2))),(df[j]/2))/
          ((1+lmd[j])*beta0[j]);

        }

      }
      if(std::strncmp(denscaller[1],"d",1)==0)
      {
        // NOTE: This is the gradient with respect to the log density
        NumericVector I0(n1),I(n1), Sign(n1),beta0(n1);

        for(int a=0;a<n1;a++){
          I0[a]=(y[a]<=mu[a]);
          I[a]= 1-I0[a];
          Sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];}
        for(j=0;j<n1;j++)
        {

          outd[j] = -(1+df[j])*(mu[j]-y[j])/(pow((mu[j]-y[j]),2)+pow(phi[j],2)*df[j]*Sign[j]);

        }


      }
    }


    else if(parCaller == "df")
    {

      if(std::strncmp(denscaller[i],"u",1)==0)
      {
        double ibeta0;
        NumericVector I0(n1),I(n1), sign(n1),sign2(n1),Z(n1),beta0(n1),digamma0,digamma1,ghy0(n1);
        NumericMatrix A(n1,3), B(n1,2),ghy(n1,1);
        for(int a=0;a<n1;a++){
          I0[a]=(y[a]<=mu[a]);
          I[a]= 1-I0[a];
          sign[a]=I0[a]*(-1) + I[a]*1;
          sign2[a]=I0[a]*(-1) + I[a]*1;
          beta0[a]=R::beta(df[a]/2, 1/2);
          digamma0[a]=R::digamma(df[a]/2);
          digamma1[a]=R::digamma((df[a]+1)/2);

          A(a,0)=0.5;
          A(a,1)=df[a]/2;
          A(a,2)=df[a]/2;
          B(a,0)=1+df[a/2];
          B(a,1)=1+df[a/2];
          Z[a]=(df[a]*pow(phi[a],2)*pow(sign[a],2))/(pow((y[a]-mu[a]),2)+pow(sign[a],2)*df[a]*pow(phi[a],2));
        }
        ghy = ghypergeo(A, B, Z, 10);

        for(j=0;j<n1;j++)
        {
          ibeta0=ibeta(Z[j], df[j]*0.5, 0.5, FALSE, TRUE);
          ghy0[j] = ghy(j,1);
          outu[j] = (sign[j]/(2*(1+lmd[j])*pow(df[j],2)*beta0[j]))*
            (sign2[j]*4*pow(Z[j],(df[j]/2))*ghy0[j]+(df[j]*(
                -2*(y[j]-mu[j])*sqrt(1/(pow((y[j]-mu[j]),2)+pow(sign[j],2)*df[j]*pow(phi[j],2)))*pow(Z[j],(df[j]/2))-
                sign2[j]*df[j]*ibeta0*
                (log(Z[j])-digamma0[j]+digamma1[j]))));

        }
      }

      if(std::strncmp(denscaller[i],"d",1)==0)
      {
        NumericVector I0(n1),I(n1), Sign(n1), C1(n1), C0(n1), C3(n1), DigammaM(n1);
        NumericMatrix A(n1,3), B(n1,2),ghy(n1,1);
        for(int a=0;a<n1;a++){
          I0[a]=(y[a]<=mu[a]);
          I[a]= 1-I0[a];
          Sign[a]=I0[a]*(-1) + I[a]*pow(lmd[a],2);

          C1[a]  =  pow((mu[a]-y[a]),2)+pow(phi[a],2)*df[a]*Sign[a];
          C0[a]  =  pow((mu[a]-y[a]),2)-pow(phi[a],2)*Sign[a];
          C3[a]  =  log(df[a]/(df[a]+pow((mu[a]-y[a]),2)/(pow(phi[a],2)*Sign[a])));
          DigammaM[a] = R::digamma(df[a]/2) - R::digamma((1+df[a])/2);

        }
        for(j=0;j<n1;j++)
        {outd[j]=(C0[j]/C1[j] + C3[j]-DigammaM[j])/2;}
      }
    }


    else if(parCaller == "phi")
    {
      if(std::strncmp(denscaller[i],"u",1)==0)
      {
        // NOTE: This is the gradient with respect to the CDF function(not the log form) .
        NumericVector I0(n1),I(n1), sign(n1),beta0(n1);

        for(int a=0;a<n1;a++){
          I0[a]=(y[a]<=mu[a]);
          I[a]= 1-I0[a];
          sign[a]=1*I0[a]+lmd[a]*I[a];
          beta0[a]=R::beta(df[a]/2, 1/2);}

        for(j=0;j<n1;j++)
        {outu[j]=-2*sign[j]*(y[j]-mu[j])*sqrt(1/pow((y[j]-mu[j]),2)+pow(sign[j],2)*df[j]*pow(phi[j],2))*
          pow(pow(sign[j],2)*df[j]*pow(phi[j],2)/(pow((y[j]-mu[j]),2)+pow(sign[j],2)*df[j]*pow(phi[j],2)),(df[j]/2))/
          ((1+lmd[j])*phi[j]*beta0[j]);}


      }

      if(std::strncmp(denscaller[i],"d",1)==0)
      {
        // NOTE: This is the gradient with respect to the log density
        NumericVector I0(n1),I(n1), Sign(n1),beta0(n1),C1(n1),C0(n1);


        for(int a=0;a<n1;a++){
          I0[a]=(y[a]<=mu[a]);
          I[a]= 1-I0[a];
          Sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
          C1[a]= pow((mu[a]-y[a]),2)+pow(phi[a],2)*df[a]*Sign[a];
          C0[a]= pow((mu[a]-y[a]),2)-pow(phi[a],2)*Sign[a];}

        for(j=0;j<n1;j++)
        {outd[j]=df[j]*C0[j]/phi[j]/C1[j];}
      }

    }

    else if(parCaller == "lmd")
    {

      if(std::strncmp(denscaller[i],"u",1)==0)
      {
        NumericVector I0(n1),I(n1), C1(n1), A(n1), B(n1), Sign(n1);
        double ibeta0;
        for(int a=0;a<n1;a++){
          I0[a]=(y[a]<=mu[a]);
          I[a]= 1-I0[a];
          outu=mu;
          for(int i=0;i<n1;i++){
            if(I0[i]==TRUE)
            {
              A[i]    = df[i]*pow(phi[i],2)/(pow((y[i]-mu[i]),2)+df[i]*pow(phi[i],2));
              ibeta0  = ibeta(A[i], df[i]*0.5, 0.5, FALSE, TRUE);
              outu[i] = -ibeta0/pow((1+lmd[i]),2);

            }

            if(I0[i]==FALSE)
            {
              B[i] = (pow((y[i]-mu[i]),2)+df[i]*pow(phi[i],2)*pow(lmd[i],2));
              A[i] = df[i]*pow(phi[i],2)*pow(lmd[i],2)/B[i];
              ibeta0  = ibeta(A[i], df[i]*0.5, 0.5, FALSE, FALSE);

              outu[i] = -(2*(1+lmd[i])*(y[i]-mu[i])*sqrt(1/B[i])*pow(A[i],(df[i]/2))+
                ibeta0/(pow((1+lmd[i]),2)*R::beta(df[i]/2,1/2)));

            }
          }


        }



        if(std::strncmp(denscaller[i],"d",1)==0)
        {
          NumericVector I0(n1),I(n1), C1(n1) , Sign(n1);
          for(int a=0;a<n1;a++){
            I0[a]=(y[a]<=mu[a]);
            I[a]= 1-I0[a];

            C1[a]=-((1+df[a]+df[a]*lmd[a])*pow((mu[a]-y[a]),2)-pow(lmd[a],3)*pow(phi[a],2)*df[a])/
            lmd[a]/(pow((mu[a]-y[a]),2)+pow(lmd[a],2)*pow(phi[a],2)*df[a]);
            Sign[a]= 1*I0[a] + C1[a]*I[a];
          }
          for(j=0;j<n1;j++)
          {outd[j]=-1/(1+lmd[j])*Sign[j];}


        }
      }

    }


    else
    { stop("No such parameter!");}



  }



  ///////////////////////////////poisson/////////////////////////////////////
  else if(type == "poisson")
  {
    mu  = par["mu"];
    NumericVector dpois0(n1);
    for(i=0;i<n1;i++)
    {
      dpois0[i]=R::dpois(floor(y[i]),mu[i], FALSE);
    }

    if(parCaller == "mu")
    {
      for(j=0;j<n1;j++)
      {
        outu[j] = (-dpois0[j]) ;// very interesting result
        outd[j] = (y[j]/mu[j] -1);
      }

    }

    else
    { stop("No such parameter!");}

  }

  //////////////////////////////////////////////////////////////////////////
  else{
    stop("This marginal density is not implemented!");}

  //////////////////////////////////////////////////////////////////////////
  List out = List::create(_["u"] = outu,
                          _["d"] = outd) ;
  return out;

}
