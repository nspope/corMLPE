genTestData <- function(M = 20, lambda = 10, beta = 0.5, rho = 0.25, sigma = 2, tau = 1, omega = 1, seed = 101){
	## function to generate a testing set, to see if model works (NOT to see if model gives biologically meaningful anwsers)
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
	attr(out, "theta") <- theta
	attr(out, "eta") <- eta
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
