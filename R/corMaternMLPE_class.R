
# pe(1|p1 + p2) 
# covariance is ...
# pe(1|p1 * p2) = pe(1|p1 + p2) + pe(1|p1:p2)
# pe(x|p1 + p2) 
# pe(x|p1 * p2) = pe(x|p1 + p2) + pe(x|p1:p2)
# sppe(1|p1 + p2, distances, 2)
# phpe(1|p1 + p2, phylo)

#' @export
corMaternMLPE <- function(value = c(0.1, 0.1), nu = 2, form = ~1, distances = FALSE, fixed = FALSE)
{
	attr(value, "formula") <- form
	attr(value, "fixed") <- fixed
  attr(value, "distances") <- distances
  attr(value, "nu") <- nu

	class(value) <- c("corMaternMLPE", "corStruct")

	value
}

#' @export
Initialize.corMaternMLPE <- function (object, data, ...)
{
  if (attr(object, "nu") <= 0 | length(attr(object, "nu")) != 1)
    stop ("smoothness parameter 'nu' must be > 0 in \"corMaternSpatial\" object")

  # set initial parameter values
  val <- as.vector(object)
  if (length(val) == 2) 
  {
    if (!all(val > 0))
    {
      stop ("'range' and 'stddev' must be > 0 in \"corMaternSpatial\" initial value")
    }
  } else {
    stop ("initial value for \"corMaternSpatial\" parameters of wrong dimension")
  }
  val <- log(val)
  object[] <- val

	attr(object, "covariate") <- getCovariate(object, data = data)
	attr(object, "factor") <- corFactor(object)
	attr(object, "logDet") <- logDet(object)

	object
}

#' @export
Dim.corMaternMLPE <- function(object, groups, ...)
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
getCovariate.corMaternMLPE <- function(object, data, ...)
{
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

    if (class(unprocessed_distances) == "matrix")
      unprocessed_distances <- list('1' = unprocessed_distances)

    if (class(unprocessed_distances) != "list")
        stop("Distance matrices for each level of the grouping factor must be supplied as a list")

    if(is.null(names(unprocessed_distances)) || !all(names(unprocessed_distances) %in% unique_groups))
      stop("Distance matrices must be supplied in a list with names matching levels of the grouping factor")

    covariate <- list()

    for (i in unique_groups)
    {
      # Integer labels for items within a group
      ind <- groups == i
      unique_labels <- unique(as.vector(unprocessed_labels[ind,]))
      labels <- rbind(match(unprocessed_labels[ind,1], unique_labels)-1,
                      match(unprocessed_labels[ind,2], unique_labels)-1)

      if (any(labels[1,]==labels[2,]))
        stop("Self-comparisons are not allowed")
      if (any(is.na(labels)))
        stop("Missing values in labels")

      # Cholesky factor of adjacency matrix
      Lfactor <- .getCovariate_cov_cpp(labels, length(unique_labels))

      # Distances between items in a group
      distances <- as.matrix(unprocessed_distances[[i]])

      if (is.null(rownames(distances)) | is.null(colnames(distances)))
        stop("Missing column or row names in distance matrices supplied to \"corMaternMLPE\" object")
      if (!all(unique_labels %in% rownames(distances)) ||
          !all(unique_labels %in% colnames(distances)) )
        stop("Inconsistent dimensions and/or missing labels in distance matrices supplied to \"corMaternMLPE\" object.\nEnsure that row/column names of distance matrices match the labels of items being compared.")

      # order appropriately (e.g. order by unique_labels)
      distances <- distances[match(rownames(distances), unique_labels),
                             match(colnames(distances), unique_labels)]

      if (any(is.infinite(distances)) || any(is.nan(distances)) || any(distances < 0))
        stop("Negative, missing, or non-finite distances supplied to \"corMaternMLPE\" object")

      if (any(distances[lower.tri(distances)] == 0.))
        stop("Offdiagonal distances == 0 in \"corMaternMLPE\", try adding a small positive constant that is << min(distances)")

      attr(labels, "distances") <- distances
      attr(labels, "ind") <- ind
      attr(labels, "Lfactor") <- Lfactor

      covariate[[i]] <- labels
    }

    attr(covariate, "groups") <- groups

    return (covariate)
  }
}

#' @export
coef.corMaternMLPE <- function (object, unconstrained = TRUE, ...) 
{
  if (attr(object, "fixed") && unconstrained) 
    return(numeric(0))
  val <- as.vector(object)
  if (length(val) == 0)
    return(val)
  if (!unconstrained)
    val <- exp(val)
  names(val) <- c("stddev", "range") 
  val
}

#' @export
"coef<-.corMaternMLPE" <- function (object, ..., value)
{
	if (length(value) != length(object)) 
		stop("cannot change the length of the parameter of a \"corMaternMLPE\" object after initializion.")

	object[] <- value

  # force recalculation of factor and determinant
  attr(object, "factor") <- NULL
  attr(object, "factor") <- corFactor(object)
	attr(object, "logDet") <- NULL
	attr(object, "logDet") <- logDet(object)

	return(object)
}

#' @export
as.matrix.corMaternMLPE <- function(object, ...){
	corMatrix(object, full=TRUE)
}

#' @export
corFactor.corMaternMLPE <- function(object, ...){
  covariate <- getCovariate(object)
  pars <- coef(object, unconstrained=FALSE)
  lapply(covariate, function(labels)
         .factor_cov_cpp (labels, attr(labels, "Lfactor"), 
                          .Matern(attr(object, "nu"), pars["stddev"], pars["range"], attr(labels, "distances"))))
}

#' @export
corMatrix.corMaternMLPE <- function(object, full=FALSE, ...){
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
recalc.corMaternMLPE <- function(object, conLin, ...){
  covariate <- getCovariate(object)
  factor <- corFactor(object)
  for (i in unique(attr(covariate, "groups")))
  {
    conLin[["Xy"]][attr(covariate[[i]], "ind"),] <- 
      .recalc_cov_cpp(covariate[[i]],
                      factor[[i]][["Lfactor"]],
                      factor[[i]][["Mfactor"]],
                      factor[[i]][["Sfactor"]],
                      conLin[["Xy"]][attr(covariate[[i]], "ind"),,drop=FALSE])
  }
	conLin[["logLik"]] <- conLin[["logLik"]] + logLik(object)
	conLin
}

#' @export
logDet.corMaternMLPE <- function(object, ...){
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

    logdet <- sum(sapply(corFactor(object), function(factor)
              sum(log(diag(factor[["Mfactor"]]))) - sum(log(factor[["Sfactor"]]))))

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
