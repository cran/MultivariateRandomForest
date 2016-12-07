// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// splitt
List splitt(NumericMatrix X, NumericMatrix Y, int m_feature, NumericVector Index, NumericMatrix Inv_Cov_Y, int Command, NumericVector ff);
RcppExport SEXP MultivariateRandomForest_splitt(SEXP XSEXP, SEXP YSEXP, SEXP m_featureSEXP, SEXP IndexSEXP, SEXP Inv_Cov_YSEXP, SEXP CommandSEXP, SEXP ffSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type m_feature(m_featureSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Inv_Cov_Y(Inv_Cov_YSEXP);
    Rcpp::traits::input_parameter< int >::type Command(CommandSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ff(ffSEXP);
    __result = Rcpp::wrap(splitt(X, Y, m_feature, Index, Inv_Cov_Y, Command, ff));
    return __result;
END_RCPP
}
// Node_cost
double Node_cost(NumericMatrix y, NumericMatrix Inv_Cov_Y, int Command);
RcppExport SEXP MultivariateRandomForest_Node_cost(SEXP ySEXP, SEXP Inv_Cov_YSEXP, SEXP CommandSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Inv_Cov_Y(Inv_Cov_YSEXP);
    Rcpp::traits::input_parameter< int >::type Command(CommandSEXP);
    __result = Rcpp::wrap(Node_cost(y, Inv_Cov_Y, Command));
    return __result;
END_RCPP
}
