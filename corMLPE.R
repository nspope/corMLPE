############################
## corMLPE class for nlme ##
############################
## Nate Pope -- npope@coa.edu
## version Sept. 7 2014

#### EXAMPLE USE ##
#
#> source('corMLPE.R')
#>
#> test.data <- genTestData(M = 16, lambda = 15) ## random dataset
#>
#> ## need to convert population lables into numbers, being careful to be consistant with both population label vectors
#> unique_pops <- unique( c(as.character(test.data$pop1), as.character(test.data$pop2)) )
#> test.data$pop1 <- match(as.character(test.data$pop1), unique_pops)
#> test.data$pop2 <- match(as.character(test.data$pop2), unique_pops)
#>
#> ## now fit model
#> test.fit <- lme(y ~ distz, random = ~ distz|spm, correlation = corMLPE(form = ~ pop1 + pop2 | spm), data = test.data)
#> plot(augPred(test.fit, primary = ~ distz))
#>
#> ## without random effect
#> test.fit.gls <- gls(y ~ distz, correlation = corMLPE(form = ~ pop1 + pop2), data = test.data)
#>
#> ## some of the computational load is in the setup of the correlation matrix. 
#> ## so if fitting multiple models to the same set of populations, you can re-use the structure of the correlation matrix from an earlier fit.
#>
#> Z <- corZ(test.fit.gls) #extract correlation structure from previous fit
#>
#> ## re-using correlation structure
#> system.time(gls(y ~ distz, correlation = corMLPE(value = 0.2, form = ~ pop1 + pop2, Z = Z), data = test.data))
#>   ## takes ~129 seconds on Asus Zenbook
#>
#> ## NOT re-using correlation structure
#> system.time(gls(y ~ distz, correlation = corMLPE(form = ~ pop1 + pop2), data = test.data))
#>   ## takes ~153 seconds on Asus Zenbook
#
###################

library(nlme)
library(plyr)
library(Matrix)

############ NLME-CONSTRUCTOR FOR MLPE CORRELATION MATRIX

Initialize.corMLPE <-  function (object, data, ...){
	form <- formula(object)
    if (!is.null(getGroupsFormula(form))) {
        attr(object, "groups") <- getGroups(object, form, data = data)
        attr(object, "Dim") <- Dim(object, attr(object, "groups"))
    }
    else {
        attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
    }
    attr(object, "covariate") <- getCovariate(object, data = data)
    attr(object, "Z") <- corZ(object)
    return(object)
    }

getCovariate.corMLPE <- function(object, data){
	library(plyr)
	if(is.null(attr(object, "covariate"))){
		form = formula(object)
		covForm <- getCovariateFormula(form)
	    if (!is.null(getGroupsFormula(form))) {
	        grps <- getGroups(object, data = data)
	        }
	    else grps <- rep(1, nrow(data))
	    covTerms <- attr(terms(covForm), "term.labels")
	    if(length(covTerms) != 2)
			stop("must have population labels in left side of formula")
		ndat <- data.frame(grps, p1 = as.character(data[[covTerms[1]]]), p2 = as.character(data[[covTerms[2]]]) )
		out <- dlply(ndat, .(grps), function(x) {v <- alply(as.matrix(x[,c(2:3)]), 1, identity); attributes(v) <- NULL; v})
		attributes(out) <- NULL
		}
	else {
		out <- attr(object, "covariate")
		}
	return(out)
	}

corZ <- function(object, ...){
	UseMethod("corZ")
	}
	
corZ.lme <- corZ.gls <- corZ.nlme <-  function(model){
	corZ(model$modelStruct$corStruct)
	}

corZ.corMLPE <- function(object, covariate = getCovariate(object)){ ## gives the structure of the correlation matrix
	if(is.null(attr(object, "Z"))){
		pdims <- Dim(object)
		out <- list()
		for (i in 1:pdims$M) {
			combs <- combn(1:pdims$len[i], 2)
			ax <- matrix(0, nrow = pdims$len[i], ncol = pdims$len[i])
			ax[lower.tri(ax)] <- as.numeric(apply(combs, 2, function(x) any(covariate[[i]][[x[1]]] %in% covariate[[i]][[x[2]]]) ))
			ax <- ax + t(ax)
			out[[i]] <- Matrix(ax) 
			}
		}
	else {
		out <- attr(object, "Z")
		}
	return(out)
	}
	
corMatrix.corMLPE <- function(object, Z = corZ(object), corr = TRUE, ...){ ## this will be the speed bottleneck in optimisation, should rewrite in C++
	rho <- as.vector(object)
	rho <- rMLPEtrans(rho)
	out <- lapply(Z, function(z) {
		z <- z * rho
		diag(z) <- 1
		if(!corr){
			z <- solve(t(chol(z))) ## this Cholesky solve is the bottleneck. Should outsource to C++ (or covert to sparseMatrix and use Matrix package)
			}
		return(z)
		})
	if(!corr) {
		lD <- unlist(lapply(out, function(z) sum(log(diag(z))) ))
		}
	else {
		lD <- NA
		}
	if (length(out) == 1)
		out <- out[[1]]
	attr(out, "logDet") <- sum(lD)
	return(out)
	}
	
coef.corMLPE <- function (object, unconstrained = TRUE, ...) {
    if (unconstrained) {
        if (attr(object, "fixed")) {
            return(numeric(0))
        }
        else {
            return(as.vector(object))
        }
    }
    rho <- as.vector(object)
    rho <- rMLPEtrans(rho)
	names(rho) <- "Rho"
	rho
    }

corMLPE <- function(value = 0.1, form = ~1, fixed = FALSE, Z = NULL){
	if(value >= 0.5 | value <= 0)
		stop("parameter must be between 0 and 0.5")
	value <- MLPEtrans(value)
	attr(value, "formula") <- form
	attr(value, "fixed") <- fixed
	if(!is.null(Z))
		attr(value, "Z") <- Z ## should add a check for this
	class(value) <- c("corMLPE", "corStruct")
	value
	}

MLPEtrans <- function(x) log( (2 * x)/(1 - 2 * x) )
rMLPEtrans <- function(z) exp(z) / ((1 + exp(z))*2)

environment(MLPEtrans) <- asNamespace("nlme")
environment(rMLPEtrans) <- asNamespace("nlme")
environment(corMLPE) <- asNamespace("nlme")
environment(corZ) <- asNamespace("nlme")
environment(corZ.corMLPE) <- asNamespace("nlme")
environment(coef.corMLPE) <- asNamespace("nlme")
environment(corMatrix.corMLPE) <- asNamespace("nlme")
environment(getCovariate.corMLPE) <- asNamespace("nlme")

genTestData <- function(M = 20, lambda = 10, beta = 0.5, rho = 0.25, sigma = 2, tau = 1, omega = 1, seed = 101){
	## function to generate a testing set, to see if model works (NOT to see if model gives biologically meaningful anwsers)
	library(MASS)
	set.seed(seed)
	y <- list()
	covariate <- list()
	distz <- list()
	spm <- list()
	Nm <- rpois(M, lambda)
	eta <- rnorm(M, 0, tau)
	theta <- rnorm(M, 0, omega)
	for(i in 1:M){
		dd0 <- dist(cbind(rnorm(Nm[i]), rnorm(Nm[i])))
		dd <- as.vector(dd0)
		comb <- combn(1:Nm[i], 2)
		P <- diag(ncol(comb))
		for(j in 1:(ncol(comb)-1)){
			for(k in (j+1):ncol(comb)){
				P[j,k] <- P[k,j] <- any(comb[,j] %in% comb[,k]) * rho
				}
			}
		y[[i]] <- mvrnorm(1, dd*beta + eta[i] + dd*theta[i], sigma*P)
		covariate[[i]] <- data.frame(t(comb))
		spm[[i]] <- rep(i, ncol(comb))
		distz[[i]] <- dd
		}
	y <- unlist(y)
	covariate <- Reduce(rbind, covariate)
	spm <- unlist(spm)
	distz <- unlist(distz)
	out <- data.frame(y=y,pop1=paste(spm,covariate[,1],sep=":"),pop2=paste(spm,covariate[,2],sep=":"),spm=spm,distz=distz,stringsAsFactors = FALSE)
	attr(out, "beta") <- beta
	attr(out, "rho") <- rho
	attr(out, "sigma") <- sigma
	attr(out, "corMat") <- P
	attr(out, "distMat") <- dd0
	out
	}

##### ORIGINAL CLARKE-TYPE MODEL

### implementation of Clarke's MLPE model, with addition of GLS variance-covariance matrix for BLUE (Clarke used OLS)
ClarkeLL <- function(par, V, A, y){
	## REML formula given in Clarke et al. 2002
	sigma <- par[1]
	rho <- rMLPEtrans(par[2]) # transform from unconstrained to constrained [-0.5, 0.5] form
	w <- nrow(V)
	V <- V*rho
	diag(V) <- 1
	phi <- solve(t(A) %*% solve(V) %*% A) %*% t(A) %*% solve(V) %*% y
	r <- y - A %*% phi
	ll <- 0.5 * ( log(det(V)) + t(r) %*% solve(V) %*% r * sigma^-1 + log(det( t(A) %*% solve(V) %*% A )) + (w-2)*log(sigma) )
	return(ll)
	}

BLUE <- function(Sol, V, y, A, pop1, pop2){
	rho <- rMLPEtrans(Sol$par[2])
	sigma <- Sol$par[1]
	st <- sigma*rho
	se <- sigma - 2*st
	## redefine variance-covariance matrix
	V_hat <- V * rho
	diag(V_hat) <- 1
	## BLUE
	phi <- solve(t(A) %*% solve(V_hat) %*% A) %*% t(A) %*% solve(V_hat) %*% y
	## variance covariance of blue
	vc.ols <- sigma * solve(t(A) %*% A) %*% t(A) %*% V_hat %*% A %*% solve(t(A) %*% A) ## this appears to be a sandwich estimator
	vc.ml <- solve(t(A) %*% solve(sigma*V_hat) %*% A) ## this is identical to gls results
	## output
	attr(phi, "vcov.ols") <- vc.ols
	attr(phi, "vcov.ml") <- vc.ml
	return(phi)
	}
	
MLPE <- function(y, X, V, pop1, pop2){
	X[,2:ncol(X)] <- apply( as.matrix(X[,2:ncol(X)]), 2, function(x) x - mean(x) )
	A <- X 
	Sol <- optim(c(1, 0), ClarkeLL, V = V, A = A, y = y, hessian = T)
	Est <- BLUE(Sol, V, y, A, pop1, pop2)
	VC <- c(sigma_e = Sol$par[1], rho = rMLPEtrans(Sol$par[2]), sigma_t = Sol$par[1] * rMLPEtrans(Sol$par[2]) )
	out <- list(BLUE=Est, VC=VC, Sol=Sol)
	class(out) <- "MLPE"
	return(out)
	}

summary.MLPE <- function(model){
	out <- matrix(nrow = length(model$BLUE), ncol = 3)
	out[,1] <- c(model$BLUE)
	out[,2] <- c(sqrt(diag(attr(model$BLUE, "vcov.ml"))) )
	out[,3] <- c(sqrt(diag(attr(model$BLUE, "vcov.ols"))) )
	colnames(out) <- c("Estimate", "GLS.se", "OLS.se")
	rownames(out) <- c("alpha", paste0("beta",1:(length(model$BLUE)-1)) )
	cat("MLPE model fit by REML in", model$Sol$counts[1], "iterations.\nREML value", model$Sol$value, "\n")
	cat("\nVariance components\n")
	cat("\tsigma_e (residual error) :", model$VC[1], "\n\trho (correlation) :", model$VC[2], "\n\tsigma_tau (rho*sigma_e) :", model$VC[3], "\n")
	cat("\nBest Linear Unbiased Estimators\n")
	print(out)
	}
