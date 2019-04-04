// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ed0_
Rcpp::NumericVector ed0_(int yday, double hour, double lat, double lon, double tcl, double o3, double cf, double albedo);
RcppExport SEXP _ed0_ed0_(SEXP ydaySEXP, SEXP hourSEXP, SEXP latSEXP, SEXP lonSEXP, SEXP tclSEXP, SEXP o3SEXP, SEXP cfSEXP, SEXP albedoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type yday(ydaySEXP);
    Rcpp::traits::input_parameter< double >::type hour(hourSEXP);
    Rcpp::traits::input_parameter< double >::type lat(latSEXP);
    Rcpp::traits::input_parameter< double >::type lon(lonSEXP);
    Rcpp::traits::input_parameter< double >::type tcl(tclSEXP);
    Rcpp::traits::input_parameter< double >::type o3(o3SEXP);
    Rcpp::traits::input_parameter< double >::type cf(cfSEXP);
    Rcpp::traits::input_parameter< double >::type albedo(albedoSEXP);
    rcpp_result_gen = Rcpp::wrap(ed0_(yday, hour, lat, lon, tcl, o3, cf, albedo));
    return rcpp_result_gen;
END_RCPP
}
// edpar_
std::vector<double> edpar_(Rcpp::IntegerVector yday, Rcpp::NumericVector hour, Rcpp::NumericVector lat, Rcpp::NumericVector lon, Rcpp::NumericVector tcl, Rcpp::NumericVector o3, Rcpp::NumericVector cf, Rcpp::NumericVector albedo);
RcppExport SEXP _ed0_edpar_(SEXP ydaySEXP, SEXP hourSEXP, SEXP latSEXP, SEXP lonSEXP, SEXP tclSEXP, SEXP o3SEXP, SEXP cfSEXP, SEXP albedoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type yday(ydaySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type hour(hourSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lat(latSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lon(lonSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tcl(tclSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type o3(o3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type cf(cfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type albedo(albedoSEXP);
    rcpp_result_gen = Rcpp::wrap(edpar_(yday, hour, lat, lon, tcl, o3, cf, albedo));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ed0_ed0_", (DL_FUNC) &_ed0_ed0_, 8},
    {"_ed0_edpar_", (DL_FUNC) &_ed0_edpar_, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_ed0(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
