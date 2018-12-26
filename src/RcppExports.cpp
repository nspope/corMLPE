// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// factor_cov
Rcpp::List factor_cov(const arma::umat& labels, const arma::mat& L, const arma::mat& K);
RcppExport SEXP _corMLPE_factor_cov(SEXP labelsSEXP, SEXP LSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(factor_cov(labels, L, K));
    return rcpp_result_gen;
END_RCPP
}
// recalc_cov
arma::mat recalc_cov(const arma::umat& labels, const arma::mat& L, const arma::mat& M, const arma::vec& stddev, arma::mat input);
RcppExport SEXP _corMLPE_recalc_cov(SEXP labelsSEXP, SEXP LSEXP, SEXP MSEXP, SEXP stddevSEXP, SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type stddev(stddevSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(recalc_cov(labels, L, M, stddev, input));
    return rcpp_result_gen;
END_RCPP
}
// getCovariate_cov
arma::mat getCovariate_cov(const arma::umat& labels, const arma::uword nodes);
RcppExport SEXP _corMLPE_getCovariate_cov(SEXP labelsSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(getCovariate_cov(labels, nodes));
    return rcpp_result_gen;
END_RCPP
}
// recalc
arma::mat recalc(const arma::umat& labels, const arma::uword nodes, const arma::mat& P, const arma::vec& D, arma::mat input, const double rho);
RcppExport SEXP _corMLPE_recalc(SEXP labelsSEXP, SEXP nodesSEXP, SEXP PSEXP, SEXP DSEXP, SEXP inputSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(recalc(labels, nodes, P, D, input, rho));
    return rcpp_result_gen;
END_RCPP
}
// adjacency
arma::mat adjacency(const arma::umat& labels, const arma::uword nodes);
RcppExport SEXP _corMLPE_adjacency(SEXP labelsSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(adjacency(labels, nodes));
    return rcpp_result_gen;
END_RCPP
}
// linear_index
arma::uword linear_index(arma::uvec indices, arma::uword dim);
RcppExport SEXP _corMLPE_linear_index(SEXP indicesSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_index(indices, dim));
    return rcpp_result_gen;
END_RCPP
}
// adjacency_NMLPE
Rcpp::List adjacency_NMLPE(const arma::umat& label, const arma::uvec& cluster);
RcppExport SEXP _corMLPE_adjacency_NMLPE(SEXP labelSEXP, SEXP clusterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type label(labelSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type cluster(clusterSEXP);
    rcpp_result_gen = Rcpp::wrap(adjacency_NMLPE(label, cluster));
    return rcpp_result_gen;
END_RCPP
}
// change_basis_NMLPE
arma::mat change_basis_NMLPE(const arma::mat& adj_0, const arma::mat& adj_1, const arma::mat& adj_2, const arma::mat& adj_3, const arma::mat& E_0, const arma::mat& E_1, const double e_2, const double e_3, const arma::mat& Linv);
RcppExport SEXP _corMLPE_change_basis_NMLPE(SEXP adj_0SEXP, SEXP adj_1SEXP, SEXP adj_2SEXP, SEXP adj_3SEXP, SEXP E_0SEXP, SEXP E_1SEXP, SEXP e_2SEXP, SEXP e_3SEXP, SEXP LinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type adj_0(adj_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type adj_1(adj_1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type adj_2(adj_2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type adj_3(adj_3SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type E_0(E_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type E_1(E_1SEXP);
    Rcpp::traits::input_parameter< const double >::type e_2(e_2SEXP);
    Rcpp::traits::input_parameter< const double >::type e_3(e_3SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Linv(LinvSEXP);
    rcpp_result_gen = Rcpp::wrap(change_basis_NMLPE(adj_0, adj_1, adj_2, adj_3, E_0, E_1, e_2, e_3, Linv));
    return rcpp_result_gen;
END_RCPP
}
// factor_cov_NMLPE
Rcpp::List factor_cov_NMLPE(const arma::umat& labels, const arma::uvec& comparisons, const arma::mat& L, const arma::mat& K);
RcppExport SEXP _corMLPE_factor_cov_NMLPE(SEXP labelsSEXP, SEXP comparisonsSEXP, SEXP LSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type comparisons(comparisonsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(factor_cov_NMLPE(labels, comparisons, L, K));
    return rcpp_result_gen;
END_RCPP
}
// recalc_cov_NMLPE
arma::mat recalc_cov_NMLPE(const arma::umat& labels, const arma::uvec& comparisons, const arma::mat& Linv, const arma::mat& Minv, const arma::vec& stddev, arma::mat input);
RcppExport SEXP _corMLPE_recalc_cov_NMLPE(SEXP labelsSEXP, SEXP comparisonsSEXP, SEXP LinvSEXP, SEXP MinvSEXP, SEXP stddevSEXP, SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type comparisons(comparisonsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Linv(LinvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Minv(MinvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type stddev(stddevSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(recalc_cov_NMLPE(labels, comparisons, Linv, Minv, stddev, input));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_corMLPE_factor_cov", (DL_FUNC) &_corMLPE_factor_cov, 3},
    {"_corMLPE_recalc_cov", (DL_FUNC) &_corMLPE_recalc_cov, 5},
    {"_corMLPE_getCovariate_cov", (DL_FUNC) &_corMLPE_getCovariate_cov, 2},
    {"_corMLPE_recalc", (DL_FUNC) &_corMLPE_recalc, 6},
    {"_corMLPE_adjacency", (DL_FUNC) &_corMLPE_adjacency, 2},
    {"_corMLPE_linear_index", (DL_FUNC) &_corMLPE_linear_index, 2},
    {"_corMLPE_adjacency_NMLPE", (DL_FUNC) &_corMLPE_adjacency_NMLPE, 2},
    {"_corMLPE_change_basis_NMLPE", (DL_FUNC) &_corMLPE_change_basis_NMLPE, 9},
    {"_corMLPE_factor_cov_NMLPE", (DL_FUNC) &_corMLPE_factor_cov_NMLPE, 4},
    {"_corMLPE_recalc_cov_NMLPE", (DL_FUNC) &_corMLPE_recalc_cov_NMLPE, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_corMLPE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
