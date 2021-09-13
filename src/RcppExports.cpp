// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// csdpArma
arma::dvec csdpArma(int n_p, int nconstraints_p, int nblocks_p, const arma::ivec& blocktypes_p, const arma::ivec& blocksizes_p, const Rcpp::List& C_p, const Rcpp::List& A_p, const arma::dvec& b_p, const arma::cube& car, Rcpp::Function func, const int printlevel);
RcppExport SEXP _Bayesrel_csdpArma(SEXP n_pSEXP, SEXP nconstraints_pSEXP, SEXP nblocks_pSEXP, SEXP blocktypes_pSEXP, SEXP blocksizes_pSEXP, SEXP C_pSEXP, SEXP A_pSEXP, SEXP b_pSEXP, SEXP carSEXP, SEXP funcSEXP, SEXP printlevelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_p(n_pSEXP);
    Rcpp::traits::input_parameter< int >::type nconstraints_p(nconstraints_pSEXP);
    Rcpp::traits::input_parameter< int >::type nblocks_p(nblocks_pSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type blocktypes_p(blocktypes_pSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type blocksizes_p(blocksizes_pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type C_p(C_pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type A_p(A_pSEXP);
    Rcpp::traits::input_parameter< const arma::dvec& >::type b_p(b_pSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type car(carSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type func(funcSEXP);
    Rcpp::traits::input_parameter< const int >::type printlevel(printlevelSEXP);
    rcpp_result_gen = Rcpp::wrap(csdpArma(n_p, nconstraints_p, nblocks_p, blocktypes_p, blocksizes_p, C_p, A_p, b_p, car, func, printlevel));
    return rcpp_result_gen;
END_RCPP
}
// alphaArma
double alphaArma(const arma::mat& X);
RcppExport SEXP _Bayesrel_alphaArma(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(alphaArma(X));
    return rcpp_result_gen;
END_RCPP
}
// l2Arma
double l2Arma(const arma::mat& X);
RcppExport SEXP _Bayesrel_l2Arma(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(l2Arma(X));
    return rcpp_result_gen;
END_RCPP
}
// l6Arma
double l6Arma(const arma::mat& X);
RcppExport SEXP _Bayesrel_l6Arma(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(l6Arma(X));
    return rcpp_result_gen;
END_RCPP
}
// pfaArma
List pfaArma(const arma::mat& X);
RcppExport SEXP _Bayesrel_pfaArma(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(pfaArma(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Bayesrel_csdpArma", (DL_FUNC) &_Bayesrel_csdpArma, 11},
    {"_Bayesrel_alphaArma", (DL_FUNC) &_Bayesrel_alphaArma, 1},
    {"_Bayesrel_l2Arma", (DL_FUNC) &_Bayesrel_l2Arma, 1},
    {"_Bayesrel_l6Arma", (DL_FUNC) &_Bayesrel_l6Arma, 1},
    {"_Bayesrel_pfaArma", (DL_FUNC) &_Bayesrel_pfaArma, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_Bayesrel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
