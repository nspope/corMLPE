#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// all of the below could use some major cleaning up
// in particular, packing everything in a struct would simplify stuff immensely

// [[Rcpp::export(.linear_index_NMLPE_cpp)]]
arma::uword linear_index (arma::uvec indices, arma::uword dim)
{
  // given 0-based row/column indices and dimension of square matrix, calculate
  // linear index. Always sorted such that the smaller index is the row.
  indices = arma::sort(indices);
  return dim * (dim - 1)/2 - (dim - indices.at(0)) * (dim - indices.at(0) - 1)/2 + indices.at(1) - indices.at(0) - 1;
}

//struct NMLPE
//{
//  NMLPE (const arma::umat& label, const arma::uvec& cluster) :
//    obs (label.n_cols),
//    label (label),
//    cluster (cluster),
//    nodes (label.max()+1),
//    clusters (cluster.max()+1),
//    combinations (...),
//    comparisons (obs),
//  {
// //  }
//
//
//}

// [[Rcpp::export(.adjacency_NMLPE_cpp)]]
Rcpp::List adjacency_NMLPE (const arma::umat& label, 
                            const arma::uvec& cluster)
{
  // DESCRIPTION
  //   This function calculates pieces of the inverse square root that do not change
  //   with parameters. Essentially, the most interpretable parameterization of the model uses
  //   random effects design matrices that are not of full rank (e.g. redundant). It is easy 
  //   to come up with a non-redundant parameterization, but it less interpretable. So here, 
  //   we calculate the inner product between the two bases, and of the full-rank basis
  //   with itself.

  arma::uword nodes = label.max()+1,
              clusters = cluster.max()+1,
              combinations = clusters * (clusters - 1) / 2;

  if (nodes != cluster.n_elem)
    Rcpp::stop("Dimension mismatch");
  
  // figure out singleton clusters: these need to be dropped from the non-redundant basis
  arma::uvec cluster_counts = arma::zeros<arma::uvec>(clusters);
  for (arma::uword i=0; i<nodes; ++i)
    cluster_counts(cluster(i)) += 1;
  arma::uvec singletons = arma::find(cluster_counts==1),
             singles (singletons.n_elem);
  for (arma::uword i=0; i<singletons.n_elem; ++i)
  {
    arma::uvec tmp = arma::find(cluster == singletons(i)); 
    singles(i) = tmp(0);
  }
  arma::uvec keep = arma::regspace<arma::uvec>(0, combinations+nodes-1);
  singles = arma::sort(singles, "descend");
  for (auto i : singles)
    keep.shed_row(i);

  // construct adjacency matrices
  arma::mat adj_11 = arma::zeros<arma::mat>(nodes, nodes);
  arma::mat adj_12 = arma::zeros<arma::mat>(nodes, clusters);
  arma::mat adj_13 = arma::zeros<arma::mat>(nodes, combinations);
  arma::mat adj_14 = arma::zeros<arma::mat>(nodes, clusters);
  arma::mat adj_32 = arma::zeros<arma::mat>(combinations, clusters);
  arma::mat adj_33 = arma::zeros<arma::mat>(combinations, combinations);
  arma::mat adj_34 = arma::zeros<arma::mat>(combinations, clusters);

  arma::uvec outer (2);
  arma::uvec comparisons (label.n_cols, arma::fill::zeros);
  for (arma::uword i=0; i<label.n_cols; ++i)
  {
    outer.at(0) = cluster(label.at(0,i));
    outer.at(1) = cluster(label.at(1,i));
    arma::uword comparison = linear_index(outer, clusters);
    bool self_comparison = outer.at(0) == outer.at(1);
    for (arma::uword j=0; j<2; ++j)
    {
      for (arma::uword k=0; k<2; ++k)
      {
        adj_11(label.at(j,i),label.at(k,i)) += 1.;
        if (!self_comparison) 
          adj_12(label.at(j,i), outer.at(k)) += 1.;
      }
      if (!self_comparison)
        adj_13(label.at(j,i), comparison) += 1.;
      if (self_comparison) 
        adj_14(label.at(j,i), outer.at(0)) += 1.;
      if (!self_comparison) 
        adj_32(comparison, outer.at(j)) += 1.;
    }
    if (!self_comparison) 
      adj_33(comparison, comparison) += 1.;
    comparisons(i) = self_comparison ? UINT_MAX : comparison + nodes; // save so no need to recalculate down the line
  }

  // assemble both adjacency matrices. The second is left in pieces
  // so to facilitate multiplication later on.
  arma::mat A = arma::join_vert(arma::join_horiz(adj_11, adj_13), 
                                arma::join_horiz(adj_13.t(), adj_33));

  arma::mat adj_0 = arma::join_vert(adj_11, arma::trans(adj_13)),
            adj_1 = arma::join_vert(adj_12, adj_32),
            adj_2 = arma::join_vert(adj_13, adj_33),
            adj_3 = arma::join_vert(adj_14, adj_34);

  // calculate cholesky factors
  // being careful to ignore singleton rows
  arma::mat L = arma::zeros<arma::mat>(arma::size(A)),
            Linv = arma::zeros<arma::mat>(arma::size(A));
  L.submat(keep, keep) = arma::trans(arma::chol(A.submat(keep, keep)));
  Linv.submat(keep, keep) = arma::inv(arma::trimatl(L.submat(keep, keep)));

  return Rcpp::List::create(
                            Rcpp::_["comparisons"] = comparisons,
                            Rcpp::_["L"] = L,  // lower cholesky of [reduced basis]' [reduced basis]
                            Rcpp::_["Linv"] = Linv,  // lower cholesky of inverse [reduced basis]' [reduced basis]
                            Rcpp::_["adj_0"] = adj_0, 
                            Rcpp::_["adj_1"] = adj_1, 
                            Rcpp::_["adj_2"] = adj_2, 
                            Rcpp::_["adj_3"] = adj_3 
                            );  // [reduced basis]' [redundant basis] 
}

// [[Rcpp::export(.change_basis_NMLPE_cpp)]]
arma::mat change_basis_NMLPE (
                        const arma::mat& adj_0,
                        const arma::mat& adj_1,
                        const arma::mat& adj_2,
                        const arma::mat& adj_3,
                        const arma::mat& E_0,
                        const arma::mat& E_1,
                        const double e_2,
                        const double e_3,
                        const arma::mat& Linv)
{
  // Given the full (redundant) covariance model, we calculate the reduced model.

  arma::mat A = arma::join_horiz(
                arma::join_horiz(adj_0 * E_0, adj_1 * E_1),
                arma::join_horiz(adj_2 * e_2, adj_3 * e_3));

  arma::mat K = Linv.t() * Linv * A;

  return K * K.t();
}

// [[Rcpp::export(.factor_cov_NMLPE_cpp)]]
Rcpp::List factor_cov_NMLPE (const arma::umat& labels, 
                       const arma::uvec& comparisons,
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
  //
  // this is very messy at the moment, what with the "remapping" of redundant effects.
  // it would be a lot nicer to package as a struct.

  arma::uword n_obs = labels.n_cols; // number of observations

  if (comparisons.n_elem != n_obs)
    Rcpp::stop("Dimension mismatch");

  // calculate standard deviations, sqrt(1 + Z.row(i) * K * Z.row(i).t())
  arma::vec stddev = arma::ones<arma::vec>(n_obs);
  for (arma::uword i=0; i<n_obs; ++i)
    for (arma::uword j=0; j<2; ++j)
    {
      if (comparisons.at(i) < UINT_MAX) // not a self comparison
      {
        stddev.at(i) += 2.0 * K(labels.at(j,i), comparisons.at(i));
        stddev.at(i) += 0.5 * K(comparisons.at(i), comparisons.at(i));
      }
      for (arma::uword k=0; k<2; ++k)
        stddev.at(i) += K(labels.at(j,i), labels.at(k,i));
    }
  stddev = arma::sqrt(stddev);

  // calculate central matrix of square root factorization
  arma::uvec keep = arma::find(K.diag() > 0.);
  arma::mat M = arma::zeros<arma::mat>(arma::size(K));
  arma::mat Minv = arma::zeros<arma::mat>(arma::size(K));
  arma::mat tmp = L.t() * K * L;
  arma::mat tmp2;
  if (!arma::chol(tmp2, arma::eye<arma::mat>(keep.n_elem, keep.n_elem) + tmp.submat(keep,keep)))
    Rcpp::stop ("Correlation matrix is not positive definite");
  M.submat(keep,keep) = tmp2;
  Minv.submat(keep,keep) = arma::inv(tmp2);

  return Rcpp::List::create(Rcpp::_["Mfactor"] = Minv.t(),
                            Rcpp::_["Sfactor"] = stddev);
}
 
 // [[Rcpp::export(.recalc_cov_NMLPE_cpp)]]
 arma::mat recalc_cov_NMLPE (const arma::umat& labels, 
                       const arma::uvec& comparisons,
                       const arma::mat&  Linv,     
                       const arma::mat&  Minv,     
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
               n_par = input.n_cols;  // number of vectors to be multiplied

   if (comparisons.n_elem != n_obs)
     Rcpp::stop("Dimension mismatch");
 
   // S^{1/2}y
   input.each_col() %= stddev;
 
   // Z'S^{1/2}y
   arma::mat Zty = arma::zeros<arma::mat>(Linv.n_rows, n_par);
   for (arma::uword i=0; i<n_obs; ++i)
     for (arma::uword j=0; j<n_par; ++j)
     {
       for (arma::uword k=0; k<2; ++k)
         Zty.at(labels.at(k,i),j) += input.at(i,j);
       if (comparisons.at(i) < UINT_MAX) // not a self comparison
         Zty.at(comparisons.at(i),j) += input.at(i,j);
     }
 
   arma::mat X = Linv * Zty; // L^{-1}Z'S^{1/2}y
   X -= Minv * X; // (I - M^{-1})L^{-1}Z'y
   X  = Linv.t() * X; // L^{-T}(I - M^{-1})L^{-1}Z'S^{1/2}y
 
   // S^{1/2} y - ZL^{-T}(I - M^{-1})L^{-1}Z' S^{1/2} y \equiv Cov^{-1/2} y
   for (arma::uword i=0; i<n_obs; ++i)
     for (arma::uword j=0; j<n_par; ++j)
     {
       for (arma::uword k=0; k<2; ++k)
         input.at(i,j) -= X.at(labels.at(k,i), j);
       if (comparisons.at(i) < UINT_MAX)
         input.at(i,j) -= X.at(comparisons.at(i),j);
     }
 
   return input;
 }
