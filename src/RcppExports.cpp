// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dsplitt
NumericVector dsplitt(NumericVector x, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd, LogicalVector log0);
RcppExport SEXP splittV1_2_dsplitt(SEXP xSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP, SEXP log0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type log0(log0SEXP);
    __result = Rcpp::wrap(dsplitt(x, mu, df, phi, lmd, log0));
    return __result;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP splittV1_2_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(rcpp_hello());
    return __result;
END_RCPP
}
// ibeta
double ibeta(double x, double a, double b, bool log0, bool reg);
RcppExport SEXP splittV1_2_ibeta(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP log0SEXP, SEXP regSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type log0(log0SEXP);
    Rcpp::traits::input_parameter< bool >::type reg(regSEXP);
    __result = Rcpp::wrap(ibeta(x, a, b, log0, reg));
    return __result;
END_RCPP
}
// splitt_kurtosis
NumericVector splitt_kurtosis(NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP splittV1_2_splitt_kurtosis(SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    __result = Rcpp::wrap(splitt_kurtosis(df, phi, lmd));
    return __result;
END_RCPP
}
// splitt_mean
NumericVector splitt_mean(NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP splittV1_2_splitt_mean(SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    __result = Rcpp::wrap(splitt_mean(mu, df, phi, lmd));
    return __result;
END_RCPP
}
// psplitt
NumericVector psplitt(NumericVector q, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP splittV1_2_psplitt(SEXP qSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    __result = Rcpp::wrap(psplitt(q, mu, df, phi, lmd));
    return __result;
END_RCPP
}
// qsplitt
NumericVector qsplitt(NumericVector p, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP splittV1_2_qsplitt(SEXP pSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    __result = Rcpp::wrap(qsplitt(p, mu, df, phi, lmd));
    return __result;
END_RCPP
}
// rsplitt
NumericVector rsplitt(int n, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP splittV1_2_rsplitt(SEXP nSEXP, SEXP muSEXP, SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    __result = Rcpp::wrap(rsplitt(n, mu, df, phi, lmd));
    return __result;
END_RCPP
}
// splitt_skewness
NumericVector splitt_skewness(NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP splittV1_2_splitt_skewness(SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    __result = Rcpp::wrap(splitt_skewness(df, phi, lmd));
    return __result;
END_RCPP
}
// splitt_var
NumericVector splitt_var(NumericVector df, NumericVector phi, NumericVector lmd);
RcppExport SEXP splittV1_2_splitt_var(SEXP dfSEXP, SEXP phiSEXP, SEXP lmdSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmd(lmdSEXP);
    __result = Rcpp::wrap(splitt_var(df, phi, lmd));
    return __result;
END_RCPP
}
