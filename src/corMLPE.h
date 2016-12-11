#ifndef CORMLPE_H
#define CORMLPE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat makeBlocks(arma::vec v, arma::vec& x, arma::uword p);
arma::vec multBlocks(arma::mat blocks, arma::uword p, arma::vec x);
arma::vec MultLambda(arma::vec x, arma::vec v, arma::uword p);
arma::vec MultLambdaTrio(const arma::vec& x, const arma::vec& v);
arma::mat MultLambdaGroups(arma::mat x, arma::mat v, arma::uvec n, arma::uvec p);
arma::vec eigenVals(double tau, arma::uword n);
arma::mat eigenVecs(arma::uword n);
arma::uvec eigenCount(arma::uword n);
arma::uvec matrixCount(arma::uword n);
arma::vec matrixVals(arma::vec lambda, arma::mat e);
double likelihood(arma::vec v, arma::uvec cnt, arma::vec l, arma::uvec lcnt, arma::vec x, arma::vec y, arma::uword p, double sigma);
double priorTauSigma(double tau, double sigma);
double logit(double p, double base);
double rlogit(double x, double base);
arma::vec sampleTauSigma(arma::vec lltausig, arma::uword iter, double tune, arma::mat e, arma::uvec cnt, arma::uvec lcnt, arma::vec x, arma::vec y, arma::uword p);
arma::ivec build_index(int nr, int nc);
arma::ivec build_pointer(int nr, int nc);
arma::vec fill_values(int nr, int nc, arma::ivec ind, arma::ivec p, arma::vec R);
arma::umat all_pairs(arma::uword p);

#endif /* CORMLPE.H */
