// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bessel_layer_simulation
Rcpp::List bessel_layer_simulation(const double& x, const double& y, const double& s, const double& t, Rcpp::NumericVector& a);
RcppExport SEXP _layeredBB_bessel_layer_simulation(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(bessel_layer_simulation(x, y, s, t, a));
    return rcpp_result_gen;
END_RCPP
}
// multi_bessel_layer_simulation
Rcpp::List multi_bessel_layer_simulation(const int& dim, const arma::vec& x, const arma::vec& y, const double& s, const double& t, Rcpp::NumericVector& a);
RcppExport SEXP _layeredBB_multi_bessel_layer_simulation(SEXP dimSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(multi_bessel_layer_simulation(dim, x, y, s, t, a));
    return rcpp_result_gen;
END_RCPP
}
// Brownian_bridge_path_sampler
Rcpp::NumericMatrix Brownian_bridge_path_sampler(const double& x, const double& y, const double& s, const double& t, Rcpp::NumericVector times);
RcppExport SEXP _layeredBB_Brownian_bridge_path_sampler(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(Brownian_bridge_path_sampler(x, y, s, t, times));
    return rcpp_result_gen;
END_RCPP
}
// multi_brownian_bridge
Rcpp::NumericMatrix multi_brownian_bridge(const int& dim, const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const double& s, const double& t, Rcpp::NumericVector times);
RcppExport SEXP _layeredBB_multi_brownian_bridge(SEXP dimSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(multi_brownian_bridge(dim, x, y, s, t, times));
    return rcpp_result_gen;
END_RCPP
}
// min_sampler
Rcpp::NumericVector min_sampler(const double& x, const double& y, const double& s, const double& t, const double& low_bound, const double& up_bound);
RcppExport SEXP _layeredBB_min_sampler(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP low_boundSEXP, SEXP up_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type low_bound(low_boundSEXP);
    Rcpp::traits::input_parameter< const double& >::type up_bound(up_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(min_sampler(x, y, s, t, low_bound, up_bound));
    return rcpp_result_gen;
END_RCPP
}
// min_Bessel_bridge_sampler
double min_Bessel_bridge_sampler(const double& x, const double& y, const double& s, const double& t, const double& min, const double& tau, const double& q);
RcppExport SEXP _layeredBB_min_Bessel_bridge_sampler(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP tauSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double& >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(min_Bessel_bridge_sampler(x, y, s, t, min, tau, q));
    return rcpp_result_gen;
END_RCPP
}
// min_Bessel_bridge_path_sampler
Rcpp::NumericMatrix min_Bessel_bridge_path_sampler(const double& x, const double& y, const double& s, const double& t, const double& min, const double& tau, Rcpp::NumericVector times, const bool& keep_min);
RcppExport SEXP _layeredBB_min_Bessel_bridge_path_sampler(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP tauSEXP, SEXP timesSEXP, SEXP keep_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< const bool& >::type keep_min(keep_minSEXP);
    rcpp_result_gen = Rcpp::wrap(min_Bessel_bridge_path_sampler(x, y, s, t, min, tau, times, keep_min));
    return rcpp_result_gen;
END_RCPP
}
// max_sampler
Rcpp::NumericVector max_sampler(const double& x, const double& y, const double& s, const double& t, const double& low_bound, const double& up_bound);
RcppExport SEXP _layeredBB_max_sampler(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP low_boundSEXP, SEXP up_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type low_bound(low_boundSEXP);
    Rcpp::traits::input_parameter< const double& >::type up_bound(up_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(max_sampler(x, y, s, t, low_bound, up_bound));
    return rcpp_result_gen;
END_RCPP
}
// max_Bessel_bridge_sampler
double max_Bessel_bridge_sampler(const double& x, const double& y, const double& s, const double& t, const double& max, const double& tau, const double& q);
RcppExport SEXP _layeredBB_max_Bessel_bridge_sampler(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP maxSEXP, SEXP tauSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type max(maxSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double& >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(max_Bessel_bridge_sampler(x, y, s, t, max, tau, q));
    return rcpp_result_gen;
END_RCPP
}
// max_Bessel_bridge_path_sampler
Rcpp::NumericMatrix max_Bessel_bridge_path_sampler(const double& x, const double& y, const double& s, const double& t, const double& max, const double& tau, Rcpp::NumericVector times, const bool& keep_max);
RcppExport SEXP _layeredBB_max_Bessel_bridge_path_sampler(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP maxSEXP, SEXP tauSEXP, SEXP timesSEXP, SEXP keep_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type max(maxSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< const bool& >::type keep_max(keep_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(max_Bessel_bridge_path_sampler(x, y, s, t, max, tau, times, keep_max));
    return rcpp_result_gen;
END_RCPP
}
// sigma_bar
double sigma_bar(const double& j, const double& x, const double& y, const double& s, const double& t, const double& l, const double& v);
RcppExport SEXP _layeredBB_sigma_bar(SEXP jSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP lSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type l(lSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_bar(j, x, y, s, t, l, v));
    return rcpp_result_gen;
END_RCPP
}
// sigma
double sigma(const double& j, const double& x, const double& y, const double& s, const double& t, const double& l, const double& v);
RcppExport SEXP _layeredBB_sigma(SEXP jSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP lSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type l(lSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma(j, x, y, s, t, l, v));
    return rcpp_result_gen;
END_RCPP
}
// phi_bar
double phi_bar(const double& j, const double& x, const double& y, const double& s, const double& t, const double& l, const double& v);
RcppExport SEXP _layeredBB_phi_bar(SEXP jSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP lSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type l(lSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_bar(j, x, y, s, t, l, v));
    return rcpp_result_gen;
END_RCPP
}
// phi
double phi(const double& j, const double& x, const double& y, const double& s, const double& t, const double& l, const double& v);
RcppExport SEXP _layeredBB_phi(SEXP jSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP lSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type l(lSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(phi(j, x, y, s, t, l, v));
    return rcpp_result_gen;
END_RCPP
}
// psi
double psi(const double& j, const double& x, const double& y, const double& s, const double& t, const double& min, const double& v);
RcppExport SEXP _layeredBB_psi(SEXP jSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(psi(j, x, y, s, t, min, v));
    return rcpp_result_gen;
END_RCPP
}
// chi
double chi(const double& j, const double& x, const double& y, const double& s, const double& t, const double& min, const double& v);
RcppExport SEXP _layeredBB_chi(SEXP jSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(chi(j, x, y, s, t, min, v));
    return rcpp_result_gen;
END_RCPP
}
// calc_SgammaK_intervals
Rcpp::NumericVector calc_SgammaK_intervals(const int& k, const double& x, const double& y, const double& s, const double& t, const double& l, const double& v);
RcppExport SEXP _layeredBB_calc_SgammaK_intervals(SEXP kSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP lSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type l(lSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_SgammaK_intervals(k, x, y, s, t, l, v));
    return rcpp_result_gen;
END_RCPP
}
// calc_SdeltaK_1_intervals
Rcpp::NumericVector calc_SdeltaK_1_intervals(const int& k, const double& x, const double& y, const double& s, const double& t, const double& min, const double& v);
RcppExport SEXP _layeredBB_calc_SdeltaK_1_intervals(SEXP kSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_SdeltaK_1_intervals(k, x, y, s, t, min, v));
    return rcpp_result_gen;
END_RCPP
}
// calc_SdeltaK_2_intervals
Rcpp::NumericVector calc_SdeltaK_2_intervals(const int& k, const double& x, const double& y, const double& s, const double& t, const double& min, const double& v);
RcppExport SEXP _layeredBB_calc_SdeltaK_2_intervals(SEXP kSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_SdeltaK_2_intervals(k, x, y, s, t, min, v));
    return rcpp_result_gen;
END_RCPP
}
// calc_SdeltaK_intervals
Rcpp::NumericVector calc_SdeltaK_intervals(const int& k, const double& x, const double& y, const double& s, const double& t, const double& min, const double& v);
RcppExport SEXP _layeredBB_calc_SdeltaK_intervals(SEXP kSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_SdeltaK_intervals(k, x, y, s, t, min, v));
    return rcpp_result_gen;
END_RCPP
}
// product_vector_elements
double product_vector_elements(const Rcpp::NumericVector& vect);
RcppExport SEXP _layeredBB_product_vector_elements(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(product_vector_elements(vect));
    return rcpp_result_gen;
END_RCPP
}
// gamma_coin
bool gamma_coin(const double& x, const double& y, const double& s, const double& t, const double& l, const double& v, int k);
RcppExport SEXP _layeredBB_gamma_coin(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP lSEXP, SEXP vSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type l(lSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(gamma_coin(x, y, s, t, l, v, k));
    return rcpp_result_gen;
END_RCPP
}
// gamma_coin_intervals
bool gamma_coin_intervals(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::NumericVector& s, const Rcpp::NumericVector& t, const double& l, const double& v, int k);
RcppExport SEXP _layeredBB_gamma_coin_intervals(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP lSEXP, SEXP vSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type l(lSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(gamma_coin_intervals(x, y, s, t, l, v, k));
    return rcpp_result_gen;
END_RCPP
}
// delta_coin
bool delta_coin(const double& x, const double& y, const double& s, const double& t, const double& min, const double& v, int k);
RcppExport SEXP _layeredBB_delta_coin(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP vSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(delta_coin(x, y, s, t, min, v, k));
    return rcpp_result_gen;
END_RCPP
}
// delta_coin_intervals
bool delta_coin_intervals(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::NumericVector& s, const Rcpp::NumericVector& t, const double& min, const double& v, int k);
RcppExport SEXP _layeredBB_delta_coin_intervals(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP minSEXP, SEXP vSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(delta_coin_intervals(x, y, s, t, min, v, k));
    return rcpp_result_gen;
END_RCPP
}
// inv_gauss_sampler
double inv_gauss_sampler(const double& mu, const double& lambda);
RcppExport SEXP _layeredBB_inv_gauss_sampler(SEXP muSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_gauss_sampler(mu, lambda));
    return rcpp_result_gen;
END_RCPP
}
// find_max
double find_max(const Rcpp::NumericVector vect);
RcppExport SEXP _layeredBB_find_max(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(find_max(vect));
    return rcpp_result_gen;
END_RCPP
}
// find_min
double find_min(const Rcpp::NumericVector vect);
RcppExport SEXP _layeredBB_find_min(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(find_min(vect));
    return rcpp_result_gen;
END_RCPP
}
// layered_brownian_bridge
Rcpp::NumericMatrix layered_brownian_bridge(const double& x, const double& y, const double& s, const double& t, const Rcpp::NumericVector& a, int l, const Rcpp::NumericVector& times);
RcppExport SEXP _layeredBB_layered_brownian_bridge(SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP aSEXP, SEXP lSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(layered_brownian_bridge(x, y, s, t, a, l, times));
    return rcpp_result_gen;
END_RCPP
}
// multi_layered_brownian_bridge
Rcpp::NumericMatrix multi_layered_brownian_bridge(const int& dim, const arma::vec& x, const arma::vec& y, const double& s, const double& t, const Rcpp::List& layers, Rcpp::NumericVector times);
RcppExport SEXP _layeredBB_multi_layered_brownian_bridge(SEXP dimSEXP, SEXP xSEXP, SEXP ySEXP, SEXP sSEXP, SEXP tSEXP, SEXP layersSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type layers(layersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(multi_layered_brownian_bridge(dim, x, y, s, t, layers, times));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_layeredBB_bessel_layer_simulation", (DL_FUNC) &_layeredBB_bessel_layer_simulation, 5},
    {"_layeredBB_multi_bessel_layer_simulation", (DL_FUNC) &_layeredBB_multi_bessel_layer_simulation, 6},
    {"_layeredBB_Brownian_bridge_path_sampler", (DL_FUNC) &_layeredBB_Brownian_bridge_path_sampler, 5},
    {"_layeredBB_multi_brownian_bridge", (DL_FUNC) &_layeredBB_multi_brownian_bridge, 6},
    {"_layeredBB_min_sampler", (DL_FUNC) &_layeredBB_min_sampler, 6},
    {"_layeredBB_min_Bessel_bridge_sampler", (DL_FUNC) &_layeredBB_min_Bessel_bridge_sampler, 7},
    {"_layeredBB_min_Bessel_bridge_path_sampler", (DL_FUNC) &_layeredBB_min_Bessel_bridge_path_sampler, 8},
    {"_layeredBB_max_sampler", (DL_FUNC) &_layeredBB_max_sampler, 6},
    {"_layeredBB_max_Bessel_bridge_sampler", (DL_FUNC) &_layeredBB_max_Bessel_bridge_sampler, 7},
    {"_layeredBB_max_Bessel_bridge_path_sampler", (DL_FUNC) &_layeredBB_max_Bessel_bridge_path_sampler, 8},
    {"_layeredBB_sigma_bar", (DL_FUNC) &_layeredBB_sigma_bar, 7},
    {"_layeredBB_sigma", (DL_FUNC) &_layeredBB_sigma, 7},
    {"_layeredBB_phi_bar", (DL_FUNC) &_layeredBB_phi_bar, 7},
    {"_layeredBB_phi", (DL_FUNC) &_layeredBB_phi, 7},
    {"_layeredBB_psi", (DL_FUNC) &_layeredBB_psi, 7},
    {"_layeredBB_chi", (DL_FUNC) &_layeredBB_chi, 7},
    {"_layeredBB_calc_SgammaK_intervals", (DL_FUNC) &_layeredBB_calc_SgammaK_intervals, 7},
    {"_layeredBB_calc_SdeltaK_1_intervals", (DL_FUNC) &_layeredBB_calc_SdeltaK_1_intervals, 7},
    {"_layeredBB_calc_SdeltaK_2_intervals", (DL_FUNC) &_layeredBB_calc_SdeltaK_2_intervals, 7},
    {"_layeredBB_calc_SdeltaK_intervals", (DL_FUNC) &_layeredBB_calc_SdeltaK_intervals, 7},
    {"_layeredBB_product_vector_elements", (DL_FUNC) &_layeredBB_product_vector_elements, 1},
    {"_layeredBB_gamma_coin", (DL_FUNC) &_layeredBB_gamma_coin, 7},
    {"_layeredBB_gamma_coin_intervals", (DL_FUNC) &_layeredBB_gamma_coin_intervals, 7},
    {"_layeredBB_delta_coin", (DL_FUNC) &_layeredBB_delta_coin, 7},
    {"_layeredBB_delta_coin_intervals", (DL_FUNC) &_layeredBB_delta_coin_intervals, 7},
    {"_layeredBB_inv_gauss_sampler", (DL_FUNC) &_layeredBB_inv_gauss_sampler, 2},
    {"_layeredBB_find_max", (DL_FUNC) &_layeredBB_find_max, 1},
    {"_layeredBB_find_min", (DL_FUNC) &_layeredBB_find_min, 1},
    {"_layeredBB_layered_brownian_bridge", (DL_FUNC) &_layeredBB_layered_brownian_bridge, 7},
    {"_layeredBB_multi_layered_brownian_bridge", (DL_FUNC) &_layeredBB_multi_layered_brownian_bridge, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_layeredBB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
