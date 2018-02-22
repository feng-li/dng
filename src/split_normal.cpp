#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List AsymNormPDF(NumericVector y,NumericVector mu, NumericVector sigma, NumericVector lmd, bool cum, bool sumStat)
{
  int len;
  len = y.size();
  NumericVector density(len),densMean(len),densStd(len);
  NumericVector I0(len),I(len), sign(len);

  if(cum) //cum=1 PDF
  {
    for(int a=0;a<len;a++)
    {
      I0[a]=(y[a]<=mu[a]);
      I[a]= 1-I0[a];
      sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
      density[a] = sqrt(2/3.1415926)*
        exp(-pow((y[a]-mu[a]),2)/(2*sigma[a]*sigma[a]*sign[a]))/
          ((1+lmd[a])*sigma[a]);
    }

  }

  else if(!cum) //cum=0 CDF
  {
    for(int a=0;a<len;a++)
    {
      I0[a]=(y[a]<=mu[a]);
      I[a]= 1-I0[a];
      sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
      if(y[a]<=mu[a])
      {
        density[a] =2/(1+lmd[a])*R::pnorm5(y[a],mu[a],sigma[a],1,0);

      }
      else if(y[a]>mu[a])
      {
        density[a] = 1/(1+lmd[a]) +
          2*lmd[a]/(1+lmd[a])*(R::pnorm5(y[a],mu[a],sigma[a],1,0)-1/2);
      }
    }
  }


  if(sumStat)
  {
    NumericVector PostC(len),Var(len);

    for(int a=0;a<len;a++)
    {
      PostC[a] = sqrt(2/3.1415926)*sigma[a]*(lmd[a]-1);
      densMean[a] = mu[a]+PostC[a];
      Var[a] = pow(pow(((3.1415926-2)/3.1415926*pow((lmd[a]-1),2)+lmd[a]),sigma[a]),2);
      densStd[a] = sqrt(Var[a]);
    }

  }


  List out= List::create(_["density"] = density,
                         _["densMean"] = densMean,
                         _["densStd"] = densStd ) ;

  return out;

}
