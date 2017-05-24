#' Correlation structure for symmetric relational data
#'
#' @param value Starting value for the correlation parameter
#' @param form A formula that gives indicators for each side of the pairwise comparison, and optionally a grouping factor. See 'Details'.
#' @param fixed Optional. Logical, fit model with the starting value for the correlation parameter fixed
#' @param covariate Optional, supply a covariate (see 'Details')
#' @export
corMLPE <- function(value = 0.1, form = ~1, fixed = FALSE, covariate = NULL){
	if(value >= 0.5 | value <= 0)
		stop("parameter must be between 0 and 0.5")
	value <- MLPEtrans(value)
	attr(value, "formula") <- form
	attr(value, "fixed") <- fixed
	if( !is.null(covariate) )
	    attr(value, "covariate") <- covariate
	class(value) <- c("corMLPE", "corStruct")
	value
}

#' @export
Initialize.corMLPE <-  function (object, data, ...){
	form <- formula(object)
	if (!is.null(getGroupsFormula(form))) {
		attr(object, "groups") <- getGroups(object, form, data = data)
		attr(object, "Dim") <- Dim(object, attr(object, "groups"))
	} else {
		attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
	}
	attr(object, "covariate") <- getCovariate(object, data = data)
	attr(object, "pattern") <- pattern(object)
	return(object)
}

#' @export
Dim.corMLPE <- function(object, groups, ...){
	if(missing(groups))
		return(attr(object, "Dim"))
	ugrp <- unique(groups)
	groups <- factor(groups, levels = ugrp)
	len <- table(groups)
	list(N = length(groups), M = length(len), maxLen = max(len), sunLenSq = 0, len = len, start = match(ugrp,groups) - 1L) 
}

#' @export
getCovariate.corMLPE <- function(object, data){
	if( !is.null( aux <- attr(object, "covariate") ) ){
		return( aux )
	} else {
		formulator = formula(object)
		covariateFormula <- getCovariateFormula(formulator)
		if ( !is.null(getGroupsFormula(formulator) )) {
			groups <- getGroups(object, data = data)
	        } else {
			groups <- rep(1, nrow(data))
		}
		covariateTerms <- attr(terms(covariateFormula), "term.labels")

		if(length(covariateTerms) != 2)
			stop("must include two variables that indicate the pairwise comparisons on left side of formula")

		aux <- cbind(paste0(groups, "_", as.character(data[[covariateTerms[1]]])), paste0(groups, "_", as.character(data[[covariateTerms[2]]]) ) )
		uniqueLabels <- unique(as.vector(aux)) 
		for(i in 1:2)
			aux[,i] <- match(aux[,i], uniqueLabels)
		
		aux <- t(apply(aux, 1, function(x) x[order(x)]))

		# check for duplicates
		if( anyDuplicated(aux) )
			stop("in at least one group, pairwise comparisons are duplicated")

		# TODO: check for complete set of pairwise comparisons
    for(i in unique(groups))
    {
      np <- length(unique(as.vector(aux[groups==i,])))
      np <- np*(np-1)/2
      problem <- paste0("group ", i, " has ", nrow(aux[groups==i,]), " pairwise comparisons.")
      if( nrow(aux[groups==i,]) != np )
      {
        cat(problem, "\n")
        stop("currently, data must contain a complete set of pairwise comparisons.
             If you don't need the capacity to handle large sample sizes, consider using
             the old (unsupported) script, which doesn't suffer from this problem.
             Found at:\n\tgithub.com/nspope/corMLPE_unsupported/\n")
      }
    }

		permutation <- order(groups, aux[,1], aux[,2])
		invpermutation <- (1:length(permutation))[permutation]

		aux <- list(permutation=permutation, invpermutation=invpermutation)

		
		return(aux)
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
    aux <- rMLPEtrans(aux)
    names(aux) <- "Rho"
    aux 
}

#' @export
"coef<-.corMLPE" <- function (object, ..., value){
	if (length(value) != length(object)) {
		stop("cannot change the length of the parameter of a \"corMLPE\" object")
	}
	object[] <- value
	attr(object, "spectrum") <- NULL
	attr(object, "spectrum") <- spectrum(object)
	attr(object, "logDet") <- NULL
	attr(object, "logDet") <- logDet(object)
	# update anything related to coefficient here, like Cholesky factor, log determinant
	return(object)
}

#' @export
as.matrix.corMLPE <- function(object, ...){
	corMatrix(object, full = TRUE)
}

#' @export
pattern <- function(object, ...){
	UseMethod("pattern")
}

#' @export
pattern.corMLPE <- function(object, ...){
	if( !is.null( aux <- attr(object, "pattern") ) ) {
		return(aux)
	} else {
		aux <- lapply( as.list(Dim(object)$len), function(n){
			p <- sqrt(1 + 8*n)/2 + 0.5
			if(p>3){
				lcnt <- eigenCount(p)
				cnt <- matrixCount(p)
				evec <- eigenVecs(p)
			} else {
				evec <- lcnt <- cnt <- NULL
			}
			list(lcnt=lcnt, cnt=cnt, evec=evec, p=p)
		})
		attr(aux, "elements") <- sapply(aux, function(x) x$p)
		return(aux)
	}
}

#' @export
spectrum <- function(object, ...){
	UseMethod("spectrum")
}

#' @export
spectrum.corMLPE <- function(object, ...){
	if( !is.null(aux <- attr(object, "spectrum")) ){
		return(aux)
	} else {
		rho <- rMLPEtrans(as.vector(object))
		# loop over groups 
		aux <- lapply( pattern(object), function(group){
			if(group$p > 3){
				lam <- sqrt( eigenVals(rho, group$p) )
				vals <- matrixVals(lam, group$evec)
				ld <- -sum( group$lcnt * log(lam) )
			} else {
				lam <- sqrt( c(1-rho, 1+rho*2) )
				vals <- c( NaN, 1/(3*lam[2])-1/(3*lam[1]), 1/(3*lam[2])+2/(3*lam[1]) )
				ld <- -sum( 2*log(lam[1]) + log(lam[2]) )
			}
			list(ld=ld, lam=lam, vals=vals)
		})
		attr(aux, "inverse") <- Reduce(cbind, lapply(aux, function(x) x$vals))
		attr(aux, "logDet") <- sum(sapply(aux, function(x) x$ld ))
		return(aux)
	}
}

#' @export
corFactor.corMLPE <- function(object, ...){
	## should complain!
	stop("Cholesky factor is not calculated for this corStruct class")
}

#' @export
corMatrix.corMLPE <- function(object, full=FALSE, ...){
	## should construct sparse matrix
	stop("Not currently implemented, matrix never formed explicitly")
}

#' @export
recalc.corMLPE <- function(object, conLin, ...){
	## need to order colLin by permutation, then unorder
	permutation <- getCovariate(object) # probably more efficient way to do this
	conLin[["Xy"]][permutation[[2]],] <- MultLambdaGroups(as.matrix(conLin[["Xy"]][permutation[[1]],]), attr(spectrum(object), "inverse"), Dim(object)$len, attr(pattern(object), "elements") )
	conLin[["logLik"]] <- conLin[["logLik"]] + logLik(object)
	return(conLin)
}

#' @export
logDet.corMLPE <- function(object, covariate = getCovariate(object), ...){
	if( !is.null(aux <- attr(object, "logDet")) ){
		return(aux)
	} else {
		aux <- attr( spectrum(object), "logDet" )
		return(-aux)
	}
}

MLPEtrans <- function(x) log( (2 * x)/(1 - 2 * x) )
rMLPEtrans <- function(z) exp(z) / ((1 + exp(z))*2)
