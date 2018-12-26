#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export(.factor_cov_cpp)]]
Rcpp::List factor_cov (const arma::umat& labels, 
                       const arma::mat&  L,
                       const arma::mat&  K)
{
  // DESCRIPTION
  //   Calculates necessary pieces of the square root factor ... at each update ...
  // INPUTS
  //   labels : adjacency list ...
  //   L : Cholesky factor of adjacency matrix implied by labels
  //   K : covariance matrix ...
  // RETURNS

  arma::uword n_obs = labels.n_cols, // number of observations
              n_pop = L.n_rows;      // number of random effects

  // calculate standard deviations, sqrt(1 + Z.row(i) * K * Z.row(i).t())
  arma::vec stddev = arma::ones<arma::vec>(n_obs);
  for (arma::uword i=0; i<n_obs; ++i)
    for (arma::uword k=0; k<2; ++k)
      for (arma::uword l=0; l<2; ++l)
        stddev(i) += K(labels.at(k,i),labels.at(l,i));
  stddev = arma::sqrt(stddev);

  // calculate central matrix of square root factorization
  arma::mat M; 
  if (!arma::chol(M, arma::eye<arma::mat>(n_pop, n_pop) + arma::trimatu(L.t()) * K * arma::trimatl(L)))
    Rcpp::stop ("Correlation matrix is not positive definite");

  return Rcpp::List::create(Rcpp::_["Mfactor"] = M.t(),
                            Rcpp::_["Lfactor"] = L,
                            Rcpp::_["Sfactor"] = stddev);
}


// [[Rcpp::export(.recalc_cov_cpp)]]
arma::mat recalc_cov (const arma::umat& labels, 
                      const arma::mat&  L,     
                      const arma::mat&  M,     
                      const arma::vec&  stddev,
                            arma::mat   input) 
{
  // DESCRIPTION
  //    Multiplies a matrix (on the left) by the square root of a correlation matrix
  // formed by imposing a spatial covariance on the random effects in an MLPE model.
  // The scheme is developed by considering the matrix square root defined in ...
  //
  // INPUTS
  //    labels :: transposed adjacency list, where columns are observations, rows are 
  //              two, and values are (base-0) indices of the two random intercepts
  //              per observation.
  //    L      :: lower-triangular Cholesky factor of the adjacency matrix implied by
  //              "labels". 
  //    M      :: lower-triangular Cholesky factor of the central piece of the square-
  //              root factorization.
  //    stddev :: marginal standard deviations implied by M.
  //    input  :: Vector(s) to be multiplied by the square-root.
  //
  // RETURNS
  //    Matrix of same dimensions as input, that are values of input multiplied on the
  // left by the square root of a spatial MLPE correlation matrix.

  arma::uword n_obs = labels.n_cols, // number of observations
              n_pop = L.n_rows,      // number of random effects
              n_par = input.n_cols;  // number of vectors to be multiplied

  // S^{1/2}y
  input.each_col() %= stddev;

  // Z'S^{1/2}y
  arma::mat Zty = arma::zeros<arma::mat>(n_pop, n_par);
  for (arma::uword i=0; i<n_obs; ++i)
    for (arma::uword j=0; j<n_par; ++j)
      for (arma::uword k=0; k<2; ++k)
        Zty.at(labels.at(k,i),j) += input.at(i,j);

  arma::mat X = arma::solve(arma::trimatl(L), Zty); // L^{-1}Z'S^{1/2}y

  X -= arma::solve(arma::trimatl(M), X); // (I - M^{-1})L^{-1}Z'y
  X  = arma::solve(arma::trimatu(L.t()), X); // L^{-T}(I - M^{-1})L^{-1}Z'S^{1/2}y

  // S^{1/2} y - ZL^{-T}(I - M^{-1})L^{-1}Z' S^{1/2} y \equiv Cov^{-1/2} y
  for (arma::uword i=0; i<n_obs; ++i)
    for (arma::uword j=0; j<n_par; ++j)
      for (arma::uword k=0; k<2; ++k)
        input.at(i,j) -= X.at(labels.at(k,i), j);

  return input;
}

//struct corCMLPE
//{
//  // dims
//  const uword n_groups,
//              n_observ,
//              n_labels;
//
//  const uvec  group_sizes,
//              groups;
//
//  const umat  labels;
//
//
//
//}

// [[Rcpp::export(.getCovariate_cov_cpp)]]
arma::mat getCovariate_cov (const arma::umat& labels, 
                            const arma::uword nodes)
{
  arma::mat adj = arma::zeros<arma::mat>(nodes, nodes);

  for (arma::uword i=0; i<labels.n_cols; ++i)
    adj.at(labels.at(0,i),labels.at(1,i)) += 1.;

  adj        += adj.t();
  adj.diag() += arma::sum(adj, 1);

  return arma::trans(arma::chol(adj));
}
