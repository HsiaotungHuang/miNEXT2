// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Hypergeometric
double Hypergeometric(int K, int k, int N, int n);
RcppExport SEXP _miNEXT2_Hypergeometric(SEXP KSEXP, SEXP kSEXP, SEXP NSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(Hypergeometric(K, k, N, n));
    return rcpp_result_gen;
END_RCPP
}
// fk
double fk(int k1, int k2, int m1, int m2, NumericVector x1, NumericVector y1);
RcppExport SEXP _miNEXT2_fk(SEXP k1SEXP, SEXP k2SEXP, SEXP m1SEXP, SEXP m2SEXP, SEXP x1SEXP, SEXP y1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< int >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< int >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y1(y1SEXP);
    rcpp_result_gen = Rcpp::wrap(fk(k1, k2, m1, m2, x1, y1));
    return rcpp_result_gen;
END_RCPP
}
// D_share
double D_share(NumericVector xi, NumericVector yi, double m1, double m2, double q);
RcppExport SEXP _miNEXT2_D_share(SEXP xiSEXP, SEXP yiSEXP, SEXP m1SEXP, SEXP m2SEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< double >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< double >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(D_share(xi, yi, m1, m2, q));
    return rcpp_result_gen;
END_RCPP
}
// h0_hat_cpp
double h0_hat_cpp(NumericVector pi1, NumericVector pi2, int m1, int m2s, int n1, int n2);
RcppExport SEXP _miNEXT2_h0_hat_cpp(SEXP pi1SEXP, SEXP pi2SEXP, SEXP m1SEXP, SEXP m2sSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< int >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< int >::type m2s(m2sSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(h0_hat_cpp(pi1, pi2, m1, m2s, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// h1_hat_cpp
double h1_hat_cpp(NumericVector pi1, NumericVector pi2, NumericVector xi1, NumericVector xi2, int m1, int m2, int n1, int n2);
RcppExport SEXP _miNEXT2_h1_hat_cpp(SEXP pi1SEXP, SEXP pi2SEXP, SEXP xi1SEXP, SEXP xi2SEXP, SEXP m1SEXP, SEXP m2SEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi1(xi1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi2(xi2SEXP);
    Rcpp::traits::input_parameter< int >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< int >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(h1_hat_cpp(pi1, pi2, xi1, xi2, m1, m2, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// un_abun
NumericVector un_abun(NumericVector xi, int n, int m);
RcppExport SEXP _miNEXT2_un_abun(SEXP xiSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(un_abun(xi, n, m));
    return rcpp_result_gen;
END_RCPP
}
// sh_abun
NumericVector sh_abun(NumericVector xi1, NumericVector xi2, int n1, int m1, int n2, int m2);
RcppExport SEXP _miNEXT2_sh_abun(SEXP xi1SEXP, SEXP xi2SEXP, SEXP n1SEXP, SEXP m1SEXP, SEXP n2SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xi1(xi1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi2(xi2SEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(sh_abun(xi1, xi2, n1, m1, n2, m2));
    return rcpp_result_gen;
END_RCPP
}
// fk_inc
double fk_inc(int k1, int k2, int T1, int T2, int t1, int t2, NumericVector x1, NumericVector y1);
RcppExport SEXP _miNEXT2_fk_inc(SEXP k1SEXP, SEXP k2SEXP, SEXP T1SEXP, SEXP T2SEXP, SEXP t1SEXP, SEXP t2SEXP, SEXP x1SEXP, SEXP y1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< int >::type T1(T1SEXP);
    Rcpp::traits::input_parameter< int >::type T2(T2SEXP);
    Rcpp::traits::input_parameter< int >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< int >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y1(y1SEXP);
    rcpp_result_gen = Rcpp::wrap(fk_inc(k1, k2, T1, T2, t1, t2, x1, y1));
    return rcpp_result_gen;
END_RCPP
}
// D0_rare
NumericVector D0_rare(NumericVector xi, NumericVector yi, double t1, double t2);
RcppExport SEXP _miNEXT2_D0_rare(SEXP xiSEXP, SEXP yiSEXP, SEXP t1SEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< double >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type t2(t2SEXP);
    rcpp_result_gen = Rcpp::wrap(D0_rare(xi, yi, t1, t2));
    return rcpp_result_gen;
END_RCPP
}
// h0_hat_inci
double h0_hat_inci(NumericVector pi1, NumericVector pi2, int t1, int t2, int T1, int T2);
RcppExport SEXP _miNEXT2_h0_hat_inci(SEXP pi1SEXP, SEXP pi2SEXP, SEXP t1SEXP, SEXP t2SEXP, SEXP T1SEXP, SEXP T2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< int >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< int >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< int >::type T1(T1SEXP);
    Rcpp::traits::input_parameter< int >::type T2(T2SEXP);
    rcpp_result_gen = Rcpp::wrap(h0_hat_inci(pi1, pi2, t1, t2, T1, T2));
    return rcpp_result_gen;
END_RCPP
}
// h1_hat_inci
double h1_hat_inci(NumericVector pi1, NumericVector pi2, NumericVector yi1, NumericVector yi2, int t1, int t2, int T1, int T2);
RcppExport SEXP _miNEXT2_h1_hat_inci(SEXP pi1SEXP, SEXP pi2SEXP, SEXP yi1SEXP, SEXP yi2SEXP, SEXP t1SEXP, SEXP t2SEXP, SEXP T1SEXP, SEXP T2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pi1(pi1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi2(pi2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi1(yi1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi2(yi2SEXP);
    Rcpp::traits::input_parameter< int >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< int >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< int >::type T1(T1SEXP);
    Rcpp::traits::input_parameter< int >::type T2(T2SEXP);
    rcpp_result_gen = Rcpp::wrap(h1_hat_inci(pi1, pi2, yi1, yi2, t1, t2, T1, T2));
    return rcpp_result_gen;
END_RCPP
}
// un_inci
NumericVector un_inci(NumericVector yi, int T, int t);
RcppExport SEXP _miNEXT2_un_inci(SEXP yiSEXP, SEXP TSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(un_inci(yi, T, t));
    return rcpp_result_gen;
END_RCPP
}
// sh_inci
NumericVector sh_inci(NumericVector yi1, NumericVector yi2, int T1, int t1, int T2, int t2);
RcppExport SEXP _miNEXT2_sh_inci(SEXP yi1SEXP, SEXP yi2SEXP, SEXP T1SEXP, SEXP t1SEXP, SEXP T2SEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type yi1(yi1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi2(yi2SEXP);
    Rcpp::traits::input_parameter< int >::type T1(T1SEXP);
    Rcpp::traits::input_parameter< int >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< int >::type T2(T2SEXP);
    Rcpp::traits::input_parameter< int >::type t2(t2SEXP);
    rcpp_result_gen = Rcpp::wrap(sh_inci(yi1, yi2, T1, t1, T2, t2));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _miNEXT2_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_miNEXT2_Hypergeometric", (DL_FUNC) &_miNEXT2_Hypergeometric, 4},
    {"_miNEXT2_fk", (DL_FUNC) &_miNEXT2_fk, 6},
    {"_miNEXT2_D_share", (DL_FUNC) &_miNEXT2_D_share, 5},
    {"_miNEXT2_h0_hat_cpp", (DL_FUNC) &_miNEXT2_h0_hat_cpp, 6},
    {"_miNEXT2_h1_hat_cpp", (DL_FUNC) &_miNEXT2_h1_hat_cpp, 8},
    {"_miNEXT2_un_abun", (DL_FUNC) &_miNEXT2_un_abun, 3},
    {"_miNEXT2_sh_abun", (DL_FUNC) &_miNEXT2_sh_abun, 6},
    {"_miNEXT2_fk_inc", (DL_FUNC) &_miNEXT2_fk_inc, 8},
    {"_miNEXT2_D0_rare", (DL_FUNC) &_miNEXT2_D0_rare, 4},
    {"_miNEXT2_h0_hat_inci", (DL_FUNC) &_miNEXT2_h0_hat_inci, 6},
    {"_miNEXT2_h1_hat_inci", (DL_FUNC) &_miNEXT2_h1_hat_inci, 8},
    {"_miNEXT2_un_inci", (DL_FUNC) &_miNEXT2_un_inci, 3},
    {"_miNEXT2_sh_inci", (DL_FUNC) &_miNEXT2_sh_inci, 6},
    {"_miNEXT2_rcpp_hello_world", (DL_FUNC) &_miNEXT2_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_miNEXT2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}