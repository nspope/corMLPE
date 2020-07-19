#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export(.recalc_cpp)]]
arma::mat recalc (const arma::umat& labels, 
                  const arma::uword nodes, 
                  const arma::mat&  P, 
                  const arma::vec&  D, 
                        arma::mat   input, 
                  const double      rho) 
{
  arma::vec M   = arma::sqrt(1. - D/(D + 1./rho - 2.));
  arma::mat X   = P * arma::diagmat((M-1.)/D) * P.t(),
            Zty = arma::zeros<arma::mat>(nodes, input.n_cols);

  // sum input by node
  for (arma::uword i=0; i<labels.n_cols; ++i)
    for (arma::uword j=0; j<input.n_cols; ++j)
      for (arma::uword k=0; k<2; ++k)
        Zty.at(labels.at(k,i),j) += input.at(i,j);

  X *= Zty;

  // map back onto input
  for (arma::uword i=0; i<labels.n_cols; ++i)
    for (arma::uword j=0; j<input.n_cols; ++j)
      for (arma::uword k=0; k<2; ++k)
        input.at(i,j) += X.at(labels.at(k,i), j);

  input /= sqrt(1. - 2.*rho);

  return input;
}

// [[Rcpp::export(.recalc_inverse_cpp)]]
arma::mat recalc_inverse (const arma::umat& labels, 
                  const arma::uword nodes, 
                  const arma::mat&  P, 
                  const arma::vec&  D, 
                        arma::mat   input, 
                  const double      rho) 
{
  arma::vec M   = arma::sqrt(1. - D/(D + 1./rho - 2.));
  arma::mat X   = P * arma::diagmat((M-1.)/D) * P.t(),
            Zty = arma::zeros<arma::mat>(nodes, input.n_cols);

  // sum input by node
  for (arma::uword i=0; i<labels.n_cols; ++i)
    for (arma::uword j=0; j<input.n_cols; ++j)
      for (arma::uword k=0; k<2; ++k)
        Zty.at(labels.at(k,i),j) += input.at(i,j);

  X *= Zty;

  // map back onto input
  for (arma::uword i=0; i<labels.n_cols; ++i)
    for (arma::uword j=0; j<input.n_cols; ++j)
      for (arma::uword k=0; k<2; ++k)
        input.at(i,j) += X.at(labels.at(k,i), j);

  input /= sqrt(1. - 2.*rho);

  // multiply by full correlation matrix
  Zty.zeros();
  for (arma::uword i=0; i<labels.n_cols; ++i)
    for (arma::uword j=0; j<input.n_cols; ++j)
      for (arma::uword k=0; k<2; ++k)
        Zty.at(labels.at(k,i),j) += input.at(i,j);

  Zty   *= rho;
  input *= (1-2.*rho);

  for (arma::uword i=0; i<labels.n_cols; ++i)
    for (arma::uword j=0; j<input.n_cols; ++j)
      for (arma::uword k=0; k<2; ++k)
        input.at(i,j) += Zty.at(labels.at(k,i), j);

  return input;
}

// [[Rcpp::export(.getCovariate_cpp)]]
arma::mat adjacency (const arma::umat& labels, 
                     const arma::uword nodes)
{
  arma::mat adj = arma::zeros<arma::mat>(nodes, nodes);

  for (arma::uword i=0; i<labels.n_cols; ++i)
    adj.at(labels.at(0,i),labels.at(1,i)) += 1.;

  adj        += adj.t();
  adj.diag() += arma::sum(adj, 1);

  return adj;
}
