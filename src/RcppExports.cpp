// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// IRLS
List IRLS(arma::mat& X, arma::mat& Tau, arma::mat& Gamma, arma::mat& Winit, bool verbose);
RcppExport SEXP _RHLP_IRLS(SEXP XSEXP, SEXP TauSEXP, SEXP GammaSEXP, SEXP WinitSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Winit(WinitSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(IRLS(X, Tau, Gamma, Winit, verbose));
    return rcpp_result_gen;
END_RCPP
}
// multinomialLogit
List multinomialLogit(arma::mat& W, arma::mat& X, arma::mat& Y, arma::mat& Gamma);
RcppExport SEXP _RHLP_multinomialLogit(SEXP WSEXP, SEXP XSEXP, SEXP YSEXP, SEXP GammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Gamma(GammaSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomialLogit(W, X, Y, Gamma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RHLP_IRLS", (DL_FUNC) &_RHLP_IRLS, 5},
    {"_RHLP_multinomialLogit", (DL_FUNC) &_RHLP_multinomialLogit, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_RHLP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
