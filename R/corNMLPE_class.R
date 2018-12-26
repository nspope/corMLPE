#note some stuff is commented out -- this is working just with the reduced basis, for debugging purposes

#' @export
corNMLPE <- function(#value = c(0.1, 0.1), 
                     value = c(0.1, 0.1, 0.1, 0.1, 0.1), 
                     nu = 2, form = ~1, clusters = FALSE, distances = FALSE, fixed = FALSE)
{
	attr(value, "formula") <- form
	attr(value, "fixed") <- fixed
  attr(value, "distances") <- distances
  attr(value, "clusters") <- clusters
  attr(value, "nu") <- nu

	class(value) <- c("corNMLPE", "corStruct")

	value
}

#' @export
Initialize.corNMLPE <- function (object, data, ...)
{
  if (attr(object, "nu") <= 0 | length(attr(object, "nu")) != 1)
    stop ("smoothness parameter 'nu' must be > 0 in \"corNMLPE\" object")

  # set initial parameter values
  val <- as.vector(object)
  if (length(val) == 5) 
#  if (length(val) == 2) 
  {
    if (!all(val > 0))
    {
      stop ("'range' and 'stddev' must be > 0 in \"corNMLPE\" initial value")
    }
  } else {
    stop ("initial value for \"corNMLPE\" parameters of wrong dimension (5 parameters needed)")
  }
  val <- log(val)
  object[] <- val

	attr(object, "covariate") <- getCovariate(object, data = data)
	attr(object, "factor") <- corFactor(object)
	attr(object, "logDet") <- logDet(object)

	object
}

#' @export
Dim.corNMLPE <- function(object, groups, ...)
{
	if(missing(groups))
		return(attr(object, "Dim"))

	ugrp   <- unique(groups)
	groups <- factor(groups, levels = ugrp)
	len    <- table(groups)

	list(N        = length(groups), 
       M        = length(len), 
       maxLen   = max(len), 
       sumLenSq = 0, 
       len      = len, 
       start    = match(ugrp,groups) - 1L) 
}

#' @export
getCovariate.corNMLPE <- function(object, data, ...)
{
  #TODO recheck this fuction
	if(!is.null( aux <- attr(object, "covariate")))
  {
		return (aux)
	} 
  else 
  {
		formulator <- formula(object)
		if (!is.null(getGroupsFormula(formulator))) 
			groups <- getGroups(object, data = data)
    else 
			groups <- as.factor(rep(1, nrow(data)))
    unique_groups <- unique(groups)

		covariateTerms <- attr(terms(getCovariateFormula(formula(object))), "term.labels")
    types_var <- sapply(data[,covariateTerms], typeof)
    if (!(length(covariateTerms) == 2 & all(types_var == "integer" || types_var == "character")))
			stop("'form' must include two variables that are of factor/integer/character type,
            that are labels for the items being compared in each pairwise observation.")

    unprocessed_labels <- cbind(as.character(data[,covariateTerms[1]]),
                                as.character(data[,covariateTerms[2]]))
    unprocessed_distances <- attr(object, "distances")
    unprocessed_clusters <- attr(object, "clusters")

    if (class(unprocessed_distances) == "matrix")
      unprocessed_distances <- list('1' = unprocessed_distances)

    if (class(unprocessed_distances) != "list")
        stop("Distance matrices for each level of the grouping factor must be supplied as a list")

    if(is.null(names(unprocessed_distances)) || !all(names(unprocessed_distances) %in% unique_groups))
      stop("Distance matrices must be supplied in a list with names matching levels of the grouping factor")

    if (class(unprocessed_clusters) != "list")
      unprocessed_clusters <- list('1' = unprocessed_clusters)

    if(is.null(names(unprocessed_clusters)) || !all(names(unprocessed_clusters) %in% unique_groups))
      stop("Cluster ids must be supplied in a list with names matching levels of the grouping factor")

    covariate <- list()

    for (i in unique_groups)
    {
      # TODO:
      # -check list of clusters as for distances above
      # -check that there's a cluster id for each unique label
      # -convert cluster ids to numeric mod 0
      # -check that distance matrix names and cluster ids match

      # Integer labels for items within a group
      ind <- groups == i
      unique_labels <- unique(as.vector(unprocessed_labels[ind,]))
      labels <- rbind(match(unprocessed_labels[ind,1], unique_labels)-1,
                      match(unprocessed_labels[ind,2], unique_labels)-1)

      if (any(labels[1,]==labels[2,]))
        stop("Self-comparisons are not allowed")
      if (any(is.na(labels)))
        stop("Missing values in labels")

      # Clusters of items
      clusters <- unprocessed_clusters[[i]]

      if (any(is.na(clusters)))
        stop("Missing values in cluster ids")
      if (is.null(names(clusters)))
        stop("Cluster ids supplied to \"corNMLPE\" object must be a named vector")
      if (!all(unique_labels %in% names(clusters)) || !all(names(clusters) %in% unique_labels))
        stop("Missing names in cluster ids supplied to \"corNMLPE\" object")

      unique_clusters <- unique(clusters)
      clusters <- match(clusters, unique_clusters)-1

      # Adjacency matrix, among other things that can be precomputed
      covr <- .adjacency_NMLPE_cpp(labels, clusters)

      # Distances between clusters in a group
      distances <- as.matrix(unprocessed_distances[[i]])

      if (is.null(rownames(distances)) | is.null(colnames(distances)))
        stop("Missing column or row names in distance matrices supplied to \"corNMLPE\" object")
      if (!all(unique_clusters %in% rownames(distances)) ||
          !all(unique_clusters %in% colnames(distances)) )
        stop("Inconsistent dimensions and/or missing clusters in distance matrices supplied to \"corNMLPE\" object.\nEnsure that row/column names of distance matrices match the labels of clusters being compared.")

      # order appropriately (e.g. order by unique_clusters)
      distances <- distances[match(rownames(distances), unique_clusters),
                             match(colnames(distances), unique_clusters)]

      if (any(is.infinite(distances)) || any(is.nan(distances)) || any(distances < 0))
        stop("Negative, missing, or non-finite distances supplied to \"corNMLPE\" object")

      if (any(distances[lower.tri(distances)] == 0.))
        stop("Offdiagonal distances == 0 in \"corNMLPE\", try adding a small positive constant that is << min(distances)")

      attr(labels, "distances") <- distances
      attr(labels, "clusters") <- clusters
      attr(labels, "comparisons") <- covr$comparisons
      attr(labels, "ind") <- ind
      attr(labels, "cholesky") <- covr$L
      attr(labels, "cholesky_inv") <- covr$Linv
      attr(labels, "n_labels") <- length(unique_labels)
      attr(labels, "n_clusters") <- length(unique_clusters)
      attr(labels, "adj_1") <- covr$adj_0 #TODO: clean all this up!
      attr(labels, "adj_2") <- covr$adj_1
      attr(labels, "adj_3") <- covr$adj_2
      attr(labels, "adj_4") <- covr$adj_3

      covariate[[i]] <- labels
    }

    attr(covariate, "groups") <- groups

    return (covariate)
  }
}

#' @export
coef.corNMLPE <- function (object, unconstrained = TRUE, ...) 
{
  if (attr(object, "fixed") && unconstrained) 
    return(numeric(0))
  val <- as.vector(object)
  if (length(val) == 0)
    return(val)
  if (!unconstrained)
    val <- exp(val)
  names(val) <- c("stddev_1", "stddev_2", "stddev_3", "stddev_4", "range")
#  names(val) <- c("stddev_1", "stddev_2")
  val
}

#' @export
"coef<-.corNMLPE" <- function (object, ..., value)
{
	if (length(value) != length(object)) 
		stop("cannot change the length of the parameter of a \"corNMLPE\" object after initializion.")

	object[] <- value

  # force recalculation of factor and determinant
  attr(object, "factor") <- NULL
  attr(object, "factor") <- corFactor(object)
	attr(object, "logDet") <- NULL
	attr(object, "logDet") <- logDet(object)

	return(object)
}

#' @export
as.matrix.corNMLPE <- function(object, ...){
	corMatrix(object, full=TRUE)
}

#' @export
corFactor.corNMLPE <- function(object, ...){
	if(!is.null(aux <- attr(object, "factor")))
	 return (aux)

  covariate <- getCovariate(object)
  pars <- coef(object, unconstrained=FALSE)
  lapply(covariate, function(labels)
         {
           K = .change_basis_NMLPE_cpp(attr(labels, "adj_1"), attr(labels, "adj_2"), attr(labels, "adj_3"), attr(labels, "adj_4"),
#                                    diag(attr(labels, "n_labels")) * pars["stddev_1"], 
#                                    diag(attr(labels, "n_clusters")) * 0, pars["stddev_2"], 0,
                                    diag(attr(labels, "n_labels")) * pars["stddev_1"], 
                                    .Matern(attr(object, "nu"), pars["stddev_2"], pars["range"], attr(labels, "distances")),
                                    pars["stddev_3"], pars["stddev_4"], 
                                    attr(labels, "cholesky_inv"))
         .factor_cov_NMLPE_cpp (labels, attr(labels, "comparisons"), attr(labels, "cholesky"), K)
         })
}

#' @export
corMatrix.corNMLPE <- function(object, full=FALSE, ...){
    # TODO:
    stop("Not implemented")
  # Depr below
  covariate <- getCovariate(object)
  pars <- coef(object, unconstrained=FALSE)
  out <- lapply (covariate, function(labels)
          {
            Sigma <- .Matern(attr(object, "nu"), pars["stddev"], pars["range"], attr(labels, "distances"))
            Z <- sparseMatrix(i=rep(1:ncol(labels), 2), j=c(labels[1,], labels[2,])+1, x=rep(1,2*ncol(labels)))
            colnames(Z) <- colnames(Sigma)
            list(Z=Z, Sigma=Sigma)
          })
  if (full)
  {
    # TODO:
    stop("Not implemented")
  }
  return (out)
}

#' @export
recalc.corNMLPE <- function(object, conLin, ...){
  covariate <- getCovariate(object)
  factor <- corFactor(object)
  for (i in unique(attr(covariate, "groups")))
  {
    conLin[["Xy"]][attr(covariate[[i]], "ind"),] <- 
      .recalc_cov_NMLPE_cpp(covariate[[i]],
                      attr(covariate[[i]], "comparisons"),
                      attr(covariate[[i]], "cholesky_inv"),
                      factor[[i]][["Mfactor"]],
                      factor[[i]][["Sfactor"]],
                      conLin[["Xy"]][attr(covariate[[i]], "ind"),,drop=FALSE])
  }
	conLin[["logLik"]] <- conLin[["logLik"]] + logLik(object)
	conLin
}

#' @export
logDet.corNMLPE <- function(object, ...){
	if( !is.null(aux <- attr(object, "logDet")) ){
		return(aux)
	} else {

    # By the matrix determinant lemma, we have 
    # logdet(I + ZXZ') = logdet(Z'Z + Xi) + logdet(X)
    #
    # logdet(Z'Z + Xi) = logdet(LL' + (Li'(M-I)Li)^{-1})
    #                  = logdet( L (I + (M-I)^{-1}) L' ) 
    #                  = 2 logdet(L) + logdet(I + (M - I)^{-1})
    #                  = 2 logdet(L) + logdet( (M-I)(M-I)^{-1} + (M-I)^{-1})
    #                  = 2 logdet(L) - logdet(M - I) + logdet(M)
    #
    # logdet(X) = logdet(Li' (M - I) Li) 
    #           = -2 logdet(L) + logdet (M - I)
    #
    # thus,
    #   logdet(I + ZXZ') = logdet(M)
    #
    # so, we need
    #   logdet(M) - logdet(S) = logdet(M) - \sum_i \log S_i

    # Note that some of the diagonal parts of M are 0 if there are redundant clusters
    # Also note that M is NOW INVERTED!
    logdet <- sum(sapply(corFactor(object), function(factor) {
              M <- factor[["Mfactor"]]
              M <- M[diag(M) > 0, diag(M) > 0]
              -sum(log(diag(M))) - sum(log(factor[["Sfactor"]]))
              }))

    return (logdet)
	}
}

.Matern <- function(nu, sd, range, distances)
{
  out <- matrix(0, nrow(distances), ncol(distances))
  distances <- sqrt(2*nu)/range * distances[lower.tri(distances)]
  out[lower.tri(out)] <- sd^2 * 2^(1-nu)/gamma(nu) * distances^nu * besselK(distances, nu)
  out <- out + t(out)
  diag(out) <- sd^2
  out
}
