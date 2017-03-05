#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double ibeta(double x, double a, double b, bool log0, bool reg)
{
  int HCond;
  HCond = (a>0 && b>0);
  if(HCond == 0)
  {
    stop("Incomplete beta function requires a > 0 and b > 0.");
  }

  double ibeta_log,lbeta0, pbeta0,out_log, out;

  pbeta0 = R::pbeta(x,a,b, TRUE, TRUE);
  //lbeta0 = R::lbeta(a,b);
  lbeta0 = ::Rf_lbeta(a,b);
  //pbeta0 = Rcpp::stats::P2<RTYPE,NA,T>pbeta(Rcpp::VectorBase<RTYPE,NA,T>&x,a,b,TRUE,TRUE);
  ibeta_log = pbeta0 + lbeta0;

  if(reg)
  { out_log = ibeta_log - lbeta0;}
  else
  {  out_log = ibeta_log;  }

  if(log0)
  {  out = out_log;}
  else
  {  out = exp(out_log);}

  return out;

}




