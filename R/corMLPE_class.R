#' Correlation structure for symmetric relational data
#'
#' @param value Starting value for the correlation parameter
#' @param form A formula including two variables that give the unordered pair of elements associated with each observation, and optionally a grouping factor that indicates the set to which the elements belong. See 'Details'.
#' @param fixed Optional. Logical, fit model with the starting value for the correlation parameter fixed
#' @details "Maximum likelihood population effects" (MLPE) is a correlation structure for dyadic, symmetric relational data: where each observation is a measurement for an unordered pair of elements from a set. For two (different) elements $i,j$, let \eqn{E[y_{i,j}]} be the expectation of the response variable (perhaps conditional on some random effects), and \deqn{y_{i,j} = E[y_{i,j}] + \alpha_{i} + \alpha_{j} + \epsilon_{i,j},} where \eqn{\alpha} are associated with unique elements of the set and are i.i.d zero-mean Gaussian random variables with standard deviation \eqn{\tau}; and \eqn{\epsilon} are i.i.d Gaussian errors with standard deviation \eqn{\sigma}. Marginally (after integrating out \eqn{\alpha}), the covariance between two observations \eqn{y_{i,j}} and \eqn{y_{k,l}} is \deqn{cov(y_{i,j}, y_{k,l}) = \tau^2 (\delta(i,k) + \delta(j,l))} where the function \eqn{\delta} evaluates to 1 when its arguments are equal and zero otherwise, and we order the indices so that \eqn{i < j, k < l} for convenience. 
#' 
#' The marginal variance is \eqn{var(y_{i,j}) = 2\tau^2 + \sigma^2}. The corresponding correlation structure has a single parameter, \eqn{\rho = \tau^2 / (2\tau^2 + \sigma^2)} which is constrained to lie between 0 and 0.5.
#'
#' The "form" argument of a corMLPE object must contain two variables that indicate the pair of elements associated with each observation, and can optionally contain a grouping factor that indicates the set to which the elements belong. Elements from different sets are treated as distinct even if they have the same label, and thus there is always a zero correlation between measurements across different sets. 
#'
#' For example, if "data.frame(elem1 = c(1,1,2), elem2 = c(2,3,4), grp = c(1,1,2))" were used as data with "form=~elem1 + elem2 | grp", then the first two observations would be correlated (because they are from the same group and share the element "1"), but would both be uncorrelated with the third observation (as the third observation is associated with the second set, despite involving an element "2" that is labelled identically to an element from the first set). The ordering within a pair does not matter. Multiple observations of the same pair of elements are allowed, as are missing combinations of pairs, but "self" comparisons are not (where both elements of a pair are the same).
#'
#' It is important to note that this correlation structure does not directly incorporate a (dis)similarity metric (which could instead be included as a covariate in the regression model), but instead tries to account for the dependence between pairwise measurements taken between the same objects.
#'
#' @references 
#' Clarke et al. 2002. Confidence limits for regression relationships between distance matrices: estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7: 361-372.
#' @examples #placeholder ... TODO
#' @export
corMLPE <- function(value = 0.1, form = ~1, fixed = FALSE){
	if(value >= 0.5 | value <= 0)
		stop("parameter must be between 0 and 0.5")

	value                  <- .MLPEtrans(value)
	attr(value, "formula") <- form
	attr(value, "fixed")   <- fixed
	class(value)           <- c("corMLPE", "corStruct")

	value
}

#' @export
Initialize.corMLPE <- function (object, data, ...){
	form <- formula(object)

	if (!is.null(getGroupsFormula(form))) 
  {
		attr(object, "groups")  <- getGroups(object, form, data = data)
		attr(object, "Dim")     <- Dim(object, attr(object, "groups"))
	} 
  else 
  {
		attr(object, "Dim")     <- Dim(object, as.factor(rep(1, nrow(data))))
	}
	attr(object, "covariate") <- getCovariate(object, data = data)

	object
}

#' @export
Dim.corMLPE <- function(object, groups, ...){
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
getCovariate.corMLPE <- function(object, data){
	if(!is.null( aux <- attr(object, "covariate")))
  {
		return (aux)
	} 
  else 
  {
		formulator       <- formula(object)
		if (!is.null(getGroupsFormula(formulator))) 
			groups <- getGroups(object, data = data)
    else 
			groups <- as.factor(rep(1, nrow(data)))
#    groups         <- getGroups(object)#TODO
		covariateTerms <- attr(terms(getCovariateFormula(formula(object))), "term.labels")
		if(length(covariateTerms) != 2)
			stop("'form' must include two factors that indicate the members of the pair for each observation.")

    # unique group labels
    labels  <- cbind(paste(groups, as.character(data[,covariateTerms[1]]), sep="_"),
                     paste(groups, as.character(data[,covariateTerms[2]]), sep="_"))
    ulabels <- unique(as.vector(labels))
    labels  <- rbind(match(labels[,1], ulabels)-1,
                     match(labels[,2], ulabels)-1)

    if (any(labels[1,]==labels[2,]))
      stop("Measurements that are self-comparisons are not allowed")

    nodes        <- length(ulabels)
    observations <- ncol(labels)

    # adjacency matrix
    adj     <- eigen(.getCovariate_cpp(labels, nodes))

    return (list(adj=adj,
                 nodes=nodes,
                 observations=observations,
                 groups=groups,
                 labels=labels))
  }
}

#' @export
coef.corMLPE <- function (object, unconstrained = TRUE, ...) {
    if (unconstrained) {
        if (attr(object, "fixed")) {
            return(numeric(0))
        }
        else {
            return(as.vector(object))
        }
    } 
    aux <- as.vector(object)
    aux <- .rMLPEtrans(aux)
    names(aux) <- "rho"
    aux 
}

#' @export
"coef<-.corMLPE" <- function (object, ..., value){
	if (length(value) != length(object)) {
		stop("cannot change the length of the parameter of a \"corMLPE\" object")
	}
	object[]               <- value
  attr(object, "logDet") <- NULL
	attr(object, "logDet") <- logDet(object)
	return(object)
}

#' @export
as.matrix.corMLPE <- function(object, ...){
	corMatrix(object, full=TRUE)
}

#' @export
corFactor.corMLPE <- function(object, ...){
	stop("The square-root factor is not explicitly formed for this this corStruct class.")
}

#' @export
corMatrix.corMLPE <- function(object, full=FALSE, gammCompatible=TRUE, ...){

  # TODO: this does not return the correct matrix if self-observations are allowed.
  # In any case I'm not sure self comparisons make sense under this model, so have
  # prohibited.

  covariate  <- getCovariate(object)
  rho        <- coef(object, unconstrained=FALSE)
  i          <- as.vector(covariate[["labels"]])+1
  j          <- rep(1:covariate[["observations"]], each=2)
  Zt         <- sparseMatrix(i=i, j=j, x=1)
  corr       <- t(Zt) %*% Zt * rho
  diag(corr) <- 1
  if (gammCompatible)
    corr     <- as(corr, "dgCMatrix_corMLPE") # see class definition below

  if (full || length(unique(covariate[["groups"]])) == 1)
    return (corr)
  else
  {
    corrl <- list()
    groups <- covariate[["groups"]]
    
    # Coerce to dense matrix primitive to allow mgcv::gamm to work.
    # The downside is a massive slowdown for large datasets.
    for (gr in unique(groups))
      if (gammCompatible)
        corrl[[gr]] <- as.matrix(corr[groups==gr,groups==gr])
      else
        corrl[[gr]] <- corr[groups==gr,groups==gr]
    return (corrl)
  }
}

#' @export
recalc.corMLPE <- function(object, conLin, ...){
  covariate    <- getCovariate(object)
  rho          <- coef(object, unconstrained=FALSE)
  conLin[["Xy"]][] <- .recalc_cpp(covariate[["labels"]],
                                  covariate[["nodes"]],
                                  covariate[["adj"]][["vectors"]],
                                  covariate[["adj"]][["values"]],
                                  conLin[["Xy"]],
                                  rho)[]
	conLin[["logLik"]] <- conLin[["logLik"]] + logLik(object)
	return(conLin)
}

#' @export
logDet.corMLPE <- function(object, covariate = getCovariate(object), ...){
	if( !is.null(aux <- attr(object, "logDet")) ){
		return(aux)
	} else {
    rho <- coef(object, unconstrained=FALSE)
    out <- covariate[["nodes"]]*log(rho/(1-2*rho)) + 
               sum(log(covariate[["adj"]][["values"]] + (1-2*rho)/rho)) +
               covariate[["observations"]]*log(1-2*rho)
    return (0.5*out)
	}
}

.MLPEtrans <- function(x) log((2 * x)/(1 - 2 * x))
.rMLPEtrans <- function(z) exp(z)/((1 + exp(z))*2)


# For mgcv::gamm,
# create a new class that inherits from dgCMatrix, e.g. essentially a dummy class,
# and overload a few methods needed for the code to work with gamm.
# If no grouping factor is included, the sparse matrix is maintained throughout
# the fitting done in gamm, with all the computational benefits that implies.
# If a grouping factor is included, gamm coerces to a dense matrix primitive.
# In this latter case the models will fit but will fart out if the dataset is
# too large. Nothing to be done about it without modifying mgcv source.

# so as not to f*** up other methods that rely on is.matrix(dgCMatrix) == FALSE
setClass("dgCMatrix_corMLPE", contains="dgCMatrix")

#' @export
is.matrix.dgCMatrix_corMLPE <- function(object, ...) TRUE
#needed to allow mgcv:::formXtViX to work

#' @export
is.numeric.dgCMatrix_corMLPE <- function(object, ...) TRUE 
#needed to allow mgcv:::formXtViX to work

#' @export
setMethod("t", signature(x = "dgCMatrix_corMLPE"),
              function(x) 
                as(.Call("Csparse_transpose", x, FALSE, PACKAGE="Matrix"), "dgCMatrix_corMLPE")
              ,
              valueClass = "dgCMatrix_corMLPE")
#needed to retain class during transpose, for mgcv:::extract.lme.cov2 to work

#' @export
setMethod("[", signature(x = "dgCMatrix_corMLPE", i = "index", j = "missing",
                         drop = "logical"),
          function (x, i,j, ..., drop) {
            cld <- getClassDef(class(x))
            na <- nargs()
            x <- if(na == 4) as(x, "dgCMatrix")[i, , drop=drop]
              else if(na == 3) as(x, "dgCMatrix")[i, drop=drop]
              else ## should not happen
                stop("Matrix-internal error in <sparseM>[i,,d]; please report")
            as(x, "dgCMatrix_corMLPE")
          })
#needed to retain class during row-subset

#' @export
setMethod("[", signature(x = "dgCMatrix_corMLPE", i = "missing", j = "index",
                         drop = "logical"),
          function (x,i,j, ..., drop) {
            cld <- getClassDef(class(x))
            x <- as(x, "dgCMatrix")[, j, drop=drop]
            as(x, "dgCMatrix_corMLPE")
          })
#needed to retain class during col-subset
