// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// gradFun4delta
NumericVector gradFun4delta(NumericMatrix u, NumericVector theta, NumericVector delta);
RcppExport SEXP _dng_gradFun4delta(SEXP uSEXP, SEXP thetaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(gradFun4delta(u, theta, delta));
    return rcpp_result_gen;
END_RCPP
}
// gradFun4df
NumericVector gradFun4df(int i, NumericMatrix rho, NumericVector df, NumericMatrix u_quantile);
RcppExport SEXP _dng_gradFun4df(SEXP iSEXP, SEXP rhoSEXP, SEXP dfSEXP, SEXP u_quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u_quantile(u_quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(gradFun4df(i, rho, df, u_quantile));
    return rcpp_result_gen;
END_RCPP
}
// gradFun4theta
NumericVector gradFun4theta(NumericMatrix u, NumericVector theta, NumericVector delta);
RcppExport SEXP _dng_gradFun4theta(SEXP uSEXP, SEXP thetaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(gradFun4theta(u, theta, delta));
    return rcpp_result_gen;
END_RCPP
}
// dCpl
NumericVector dCpl(std::string CplNM, NumericMatrix u, List parCpl, bool log0);
RcppExport SEXP _dng_dCpl(SEXP CplNMSEXP, SEXP uSEXP, SEXP parCplSEXP, SEXP log0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type CplNM(CplNMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< List >::type parCpl(parCplSEXP);
    Rcpp::traits::input_parameter< bool >::type log0(log0SEXP);
    rcpp_result_gen = Rcpp::wrap(dCpl(CplNM, u, parCpl, log0));
    return rcpp_result_gen;
END_RCPP
}
// dmvNormVecFun
NumericVector dmvNormVecFun(int i, NumericVector x, NumericVector rho);
RcppExport SEXP _dng_dmvNormVecFun(SEXP iSEXP, SEXP xSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvNormVecFun(i, x, rho));
    return rcpp_result_gen;
END_RCPP
}
// dng_grad
List dng_grad(NumericVector y, List par, std::string type, std::string parCaller, GenericVector denscaller);
RcppExport SEXP _dng_dng_grad(SEXP ySEXP, SEXP parSEXP, SEXP typeSEXP, SEXP parCallerSEXP, SEXP denscallerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type par(parSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< std::string >::type parCaller(parCallerSEXP);
    Rcpp::traits::input_parameter< GenericVector >::type denscaller(denscallerSEXP);
    rcpp_result_gen = Rcpp::wrap(dng_grad(y, par, type, parCaller, denscaller));
    return rcpp_result_gen;
END_RCPP
}
// dsplitn
NumericVector dsplitn(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector lmd, bool logarithm);
RcppExport SEXP _dng_dsplitn(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lmdSEXP, SEXP logarithmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    Rcpp::traits::input_parameter< bool >::type logarithm(logarithmSEXP);
    rcpp_result_gen = Rcpp::wrap(dsplitn(x, mu, sigma, lmd, logarithm));
    return rcpp_result_gen;
END_RCPP
}
// dsplitt
NumericVector dsplitt(NumericVector x, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd, bool logarithm);
RcppExport SEXP _dng_dsplitt(SEXP xSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP, SEXP logarithmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    Rcpp::traits::input_parameter< bool >::type logarithm(logarithmSEXP);
    rcpp_result_gen = Rcpp::wrap(dsplitt(x, mu, df, phi, lmd, logarithm));
    return rcpp_result_gen;
END_RCPP
}
// FUN
NumericMatrix FUN(int i, NumericMatrix x, NumericVector mu, NumericVector df, NumericMatrix rho, NumericVector uIdx);
RcppExport SEXP _dng_FUN(SEXP iSEXP, SEXP xSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP rhoSEXP, SEXP uIdxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type uIdx(uIdxSEXP);
    rcpp_result_gen = Rcpp::wrap(FUN(i, x, mu, df, rho, uIdx));
    return rcpp_result_gen;
END_RCPP
}
// ghypergeo
NumericVector ghypergeo(NumericMatrix a, NumericMatrix b, NumericVector z, int k);
RcppExport SEXP _dng_ghypergeo(SEXP aSEXP, SEXP bSEXP, SEXP zSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(ghypergeo(a, b, z, k));
    return rcpp_result_gen;
END_RCPP
}
// ibeta
double ibeta(double x, double a, double b, bool log0, bool reg);
RcppExport SEXP _dng_ibeta(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP log0SEXP, SEXP regSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type log0(log0SEXP);
    Rcpp::traits::input_parameter< bool >::type reg(regSEXP);
    rcpp_result_gen = Rcpp::wrap(ibeta(x, a, b, log0, reg));
    return rcpp_result_gen;
END_RCPP
}
// splitt_kurtosis
NumericVector splitt_kurtosis(NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP _dng_splitt_kurtosis(SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(splitt_kurtosis(df, phi, lmd));
    return rcpp_result_gen;
END_RCPP
}
// logCplGrad
List logCplGrad(std::string CplNM, NumericMatrix u, List parCpl, std::string parCaller);
RcppExport SEXP _dng_logCplGrad(SEXP CplNMSEXP, SEXP uSEXP, SEXP parCplSEXP, SEXP parCallerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type CplNM(CplNMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< List >::type parCpl(parCplSEXP);
    Rcpp::traits::input_parameter< std::string >::type parCaller(parCallerSEXP);
    rcpp_result_gen = Rcpp::wrap(logCplGrad(CplNM, u, parCpl, parCaller));
    return rcpp_result_gen;
END_RCPP
}
// logDensFun
NumericVector logDensFun(NumericMatrix u, NumericVector theta, NumericVector delta);
RcppExport SEXP _dng_logDensFun(SEXP uSEXP, SEXP thetaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(logDensFun(u, theta, delta));
    return rcpp_result_gen;
END_RCPP
}
// splitn_mean
NumericVector splitn_mean(NumericVector mu, NumericVector sigma, NumericVector lmd);
RcppExport SEXP _dng_splitn_mean(SEXP muSEXP, SEXP sigmaSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(splitn_mean(mu, sigma, lmd));
    return rcpp_result_gen;
END_RCPP
}
// splitt_mean
NumericVector splitt_mean(NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP _dng_splitt_mean(SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(splitt_mean(mu, df, phi, lmd));
    return rcpp_result_gen;
END_RCPP
}
// pochhammer
NumericMatrix pochhammer(NumericVector a, IntegerVector n, bool log0);
RcppExport SEXP _dng_pochhammer(SEXP aSEXP, SEXP nSEXP, SEXP log0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type log0(log0SEXP);
    rcpp_result_gen = Rcpp::wrap(pochhammer(a, n, log0));
    return rcpp_result_gen;
END_RCPP
}
// psplitn
NumericVector psplitn(NumericVector q, NumericVector mu, NumericVector sigma, NumericVector lmd, bool logarithm);
RcppExport SEXP _dng_psplitn(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lmdSEXP, SEXP logarithmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    Rcpp::traits::input_parameter< bool >::type logarithm(logarithmSEXP);
    rcpp_result_gen = Rcpp::wrap(psplitn(q, mu, sigma, lmd, logarithm));
    return rcpp_result_gen;
END_RCPP
}
// psplitt
NumericVector psplitt(NumericVector q, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP _dng_psplitt(SEXP qSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(psplitt(q, mu, df, phi, lmd));
    return rcpp_result_gen;
END_RCPP
}
// qsplitn
NumericVector qsplitn(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector lmd);
RcppExport SEXP _dng_qsplitn(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(qsplitn(p, mu, sigma, lmd));
    return rcpp_result_gen;
END_RCPP
}
// qsplitt
NumericVector qsplitt(NumericVector p, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP _dng_qsplitt(SEXP pSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(qsplitt(p, mu, df, phi, lmd));
    return rcpp_result_gen;
END_RCPP
}
// rsplitn
NumericVector rsplitn(int n, NumericVector mu, NumericVector sigma, NumericVector lmd);
RcppExport SEXP _dng_rsplitn(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(rsplitn(n, mu, sigma, lmd));
    return rcpp_result_gen;
END_RCPP
}
// rsplitt
NumericVector rsplitt(int n, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP _dng_rsplitt(SEXP nSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(rsplitt(n, mu, df, phi, lmd));
    return rcpp_result_gen;
END_RCPP
}
// splitt_skewness
NumericVector splitt_skewness(NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP _dng_splitt_skewness(SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(splitt_skewness(df, phi, lmd));
    return rcpp_result_gen;
END_RCPP
}
// splitn_kurtosis
NumericVector splitn_kurtosis(NumericVector lmd);
RcppExport SEXP _dng_splitn_kurtosis(SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(splitn_kurtosis(lmd));
    return rcpp_result_gen;
END_RCPP
}
// splitn_skewness
NumericVector splitn_skewness(NumericVector sigma, NumericVector lmd);
RcppExport SEXP _dng_splitn_skewness(SEXP sigmaSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(splitn_skewness(sigma, lmd));
    return rcpp_result_gen;
END_RCPP
}
// splitn_var
NumericVector splitn_var(NumericVector sigma, NumericVector lmd);
RcppExport SEXP _dng_splitn_var(SEXP sigmaSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(splitn_var(sigma, lmd));
    return rcpp_result_gen;
END_RCPP
}
// splitt_var
NumericVector splitt_var(NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP _dng_splitt_var(SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    rcpp_result_gen = Rcpp::wrap(splitt_var(df, phi, lmd));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dng_gradFun4delta", (DL_FUNC) &_dng_gradFun4delta, 3},
    {"_dng_gradFun4df", (DL_FUNC) &_dng_gradFun4df, 4},
    {"_dng_gradFun4theta", (DL_FUNC) &_dng_gradFun4theta, 3},
    {"_dng_dCpl", (DL_FUNC) &_dng_dCpl, 4},
    {"_dng_dmvNormVecFun", (DL_FUNC) &_dng_dmvNormVecFun, 3},
    {"_dng_dng_grad", (DL_FUNC) &_dng_dng_grad, 5},
    {"_dng_dsplitn", (DL_FUNC) &_dng_dsplitn, 5},
    {"_dng_dsplitt", (DL_FUNC) &_dng_dsplitt, 6},
    {"_dng_FUN", (DL_FUNC) &_dng_FUN, 6},
    {"_dng_ghypergeo", (DL_FUNC) &_dng_ghypergeo, 4},
    {"_dng_ibeta", (DL_FUNC) &_dng_ibeta, 5},
    {"_dng_splitt_kurtosis", (DL_FUNC) &_dng_splitt_kurtosis, 3},
    {"_dng_logCplGrad", (DL_FUNC) &_dng_logCplGrad, 4},
    {"_dng_logDensFun", (DL_FUNC) &_dng_logDensFun, 3},
    {"_dng_splitn_mean", (DL_FUNC) &_dng_splitn_mean, 3},
    {"_dng_splitt_mean", (DL_FUNC) &_dng_splitt_mean, 4},
    {"_dng_pochhammer", (DL_FUNC) &_dng_pochhammer, 3},
    {"_dng_psplitn", (DL_FUNC) &_dng_psplitn, 5},
    {"_dng_psplitt", (DL_FUNC) &_dng_psplitt, 5},
    {"_dng_qsplitn", (DL_FUNC) &_dng_qsplitn, 4},
    {"_dng_qsplitt", (DL_FUNC) &_dng_qsplitt, 5},
    {"_dng_rsplitn", (DL_FUNC) &_dng_rsplitn, 4},
    {"_dng_rsplitt", (DL_FUNC) &_dng_rsplitt, 5},
    {"_dng_splitt_skewness", (DL_FUNC) &_dng_splitt_skewness, 3},
    {"_dng_splitn_kurtosis", (DL_FUNC) &_dng_splitn_kurtosis, 1},
    {"_dng_splitn_skewness", (DL_FUNC) &_dng_splitn_skewness, 2},
    {"_dng_splitn_var", (DL_FUNC) &_dng_splitn_var, 2},
    {"_dng_splitt_var", (DL_FUNC) &_dng_splitt_var, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dng(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
