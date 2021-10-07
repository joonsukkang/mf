
# solve for tau
SOL.TAU <- function(Y, Omega, A.l, A.f, C.l, C.f){
  1/(mean(((Y-tcrossprod(A.l,A.f))[Omega])^2)
     + mean((tcrossprod(C.l, C.f))[Omega])
     - mean((tcrossprod(A.l^2,A.f^2))[Omega])
  )
}


ELBO <- function(n, p, n.obs, tau, tau.l, tau.f, A.l, A.f, C.l, C.f){
  K <- length(tau.l)
  obj <-   -n.obs/2*(log(2*pi)-log(tau)+1) +
    -1/2*(sum(C.l %*% tau.l) - n*K - n*sum(log(tau.l)) - sum(log(C.l-A.l^2))) +
    -1/2*(sum(C.f %*% tau.f) - p*K - p*sum(log(tau.f)) - sum(log(C.f-A.f^2)))
  return(obj)
}

init.K <- function(Y, K){

  n <- nrow(Y);  p <- ncol(Y)
  Omega <- !is.na(Y)


  # initialize A/B/C matrices
  A.l <- matrix(rnorm(n=n*K), nrow=n, ncol=K)
  A.f <- matrix(rnorm(n=p*K), nrow=p, ncol=K)
  B.l <- crossprod(A.l)
  B.f <- crossprod(A.f)
  C.l <- A.l^2
  C.f <- A.f^2

  # initialize tau, tau.l, tau.f
  tau <- SOL.TAU(Y, Omega, A.l, A.f, C.l, C.f)
  tau.l <- 1/colMeans(C.l);   tau.f <- 1/colMeans(C.f)

  out.list <- list(A.l=A.l, B.l=B.l, C.l=C.l,
                   A.f=A.f, B.f=B.f, C.f=C.f,
                   tau=tau, tau.l=tau.l, tau.f=tau.f)
  return(out.list)
}

# update q.l (A.l and C.l)
SOLQL.R <- function(Y, A.f, C.f, Omega, tau.l, tau){

  n <- nrow(Y)
  K <- length(tau.l)
  A.l <- matrix(0, nrow=n, ncol=K)
  C.l <- matrix(0, nrow=n, ncol=K)

  for (i in 1:n){
    A.f.i <- A.f[Omega[i,],, drop=FALSE]
    C.f.i <- C.f[Omega[i,],, drop=FALSE]
    B.f.i <- crossprod(A.f.i); diag(B.f.i) <- colSums(C.f.i)

    A.l[i,] <- Y[i,Omega[i,], drop=FALSE] %*% A.f.i %*% solve(B.f.i + diag(x=tau.l/tau, nrow=K, ncol=K))
    C.l[i,] <- A.l[i,]^2 + 1/(tau.l + tau*diag(B.f.i))
  }

  out.list <- list(A.l = A.l,
                   C.l = C.l)
  return(out.list)
}



backfit <- function(Y, init.obj, maxiter=500){

  for(i in 1:length(init.obj)) assign(names(init.obj)[i], init.obj[[i]])

  n <- nrow(Y);  p <- ncol(Y)
  Omega <- !is.na(Y)
  n.obs <- sum(Omega)
  K <- length(tau.l)
  tY <- t(Y)
  tOmega <- t(Omega)


  elbo.tol <- sqrt(.Machine$double.eps) *prod(dim(Y))
  elbo.vec <- c()

  ################### update
  for (n.iter in 1:maxiter){
    # update q_l
    #obj.temp <- SOLQL.R(Y, A.f, C.f, Omega, tau.l, tau) # R version
    obj.temp <- SOLQL(Y, A.f, C.f, Omega, tau.l, tau) # RcppArmadillo version
    A.l <- obj.temp$A.l
    C.l <- obj.temp$C.l
    rm(obj.temp)

    # update q_f
    #obj.temp <- SOLQL.R(t(Y), A.l, C.l, t(Omega), tau.f, tau)
    obj.temp <- SOLQL(tY, A.l, C.l, tOmega, tau.f, tau)
    A.f <- obj.temp$A.l
    C.f <- obj.temp$C.l
    rm(obj.temp)

    # update tau, tau.l, tau.f
    #tau <- SOL.TAU(Y, Omega, A.l, A.f, C.l, C.f)
    tau <- SOL.TAU(Y, Omega, A.l, A.f, C.l, C.f)
    tau.l <- 1/colMeans(C.l)
    tau.f <- 1/colMeans(C.f)

    # compute elbo
    elbo  <- ELBO(n, p, n.obs, tau, tau.l, tau.f, A.l, A.f, C.l, C.f)
    elbo.vec <- c(elbo.vec, elbo)

    if(n.iter>1){
      if(diff(tail(elbo.vec, 2))<elbo.tol) break
    }
  }

  out.list <- list(A.l=A.l, B.l=B.l, C.l=C.l,
                   A.f=A.f, B.f=B.f, C.f=C.f,
                   tau=tau, tau.l=tau.l, tau.f=tau.f, elbo.vec=elbo.vec)
  return(out.list)
}


greedy.init <- function(Y, greedy.Kmax=50){

  n <- nrow(Y);  p <- ncol(Y)
  Omega <- !is.na(Y)
  Y.res <- Y # 'residual' Y
  n.obs <- sum(Omega)

  elbo.tol <- sqrt(.Machine$double.eps) *prod(dim(Y))
  elbo.vec <- c()

  # first factor
  fit.temp <- backfit(Y.res, init.K(Y.res, K=1), maxiter=100)
  A.l <- fit.temp$A.l; C.l <- fit.temp$C.l
  B.l <- crossprod(A.l); diag(B.l) <- colSums(C.l)

  A.f <- fit.temp$A.f; C.f <- fit.temp$C.f
  B.f <- crossprod(A.f); diag(B.f) <- colSums(C.f)

  tau <- fit.temp$tau
  tau.l <- fit.temp$tau.l; tau.f <- fit.temp$tau.f

  elbo  <- ELBO(n, p, n.obs, tau, tau.l, tau.f, A.l, A.f, C.l, C.f)
  #print(paste0(1, '; elbo=',elbo))
  elbo.vec <- c(elbo.vec, elbo)

  fitted.k <-  fit.temp$A.l%*%t(fit.temp$A.f)
  Y.res <- Y.res - fitted.k


  # add factors
  for (k in 2:greedy.Kmax){
    fit.temp <- backfit(Y.res, init.K(Y.res, K=1), maxiter=500)

    # null checking
    A.l <- cbind(A.l, fit.temp$A.l); C.l <- cbind(C.l, fit.temp$C.l)
    B.l <- crossprod(A.l); diag(B.l) <- colSums(C.l)

    A.f <- cbind(A.f, fit.temp$A.f); C.f <- cbind(C.f, fit.temp$C.f)
    B.f <- crossprod(A.f); diag(B.f) <- colSums(C.f)

    tau <- SOL.TAU(Y, Omega, A.l, A.f, C.l, C.f)
    tau.l <- c(tau.l, fit.temp$tau.l)
    tau.f <- c(tau.f, fit.temp$tau.f)

    elbo  <- ELBO(n, p, n.obs, tau, tau.l, tau.f, A.l, A.f, C.l, C.f)
    #print(paste0(k, '; elbo=',elbo))
    elbo.vec <- c(elbo.vec, elbo)

    if(diff(tail(elbo.vec,2))<elbo.tol){
      A.l <- A.l[,-k,drop=FALSE]; C.l <- C.l[,-k,drop=FALSE]; B.l <- crossprod(A.l); diag(B.l) <- colSums(C.l)
      A.f <- A.f[,-k,drop=FALSE]; C.f <- C.f[,-k,drop=FALSE]; B.f <- crossprod(A.f); diag(B.f) <- colSums(C.f)
      tau.l <- tau.l[-k]; tau.f <- tau.f[-k]
      tau <- SOL.TAU(Y, Omega, A.l, A.f, C.l, C.f)

      #print(paste0('factor ', k, ' deleted; total ', k-1, ' factors fitted'))
      break
    }

    fitted.k <-  fit.temp$A.l %*% t(fit.temp$A.f)
    Y.res <- Y.res - fitted.k
  }

  out.list <- list(A.l=A.l, B.l=B.l, C.l=C.l,
                   A.f=A.f, B.f=B.f, C.f=C.f,
                   tau=tau, tau.l=tau.l, tau.f=tau.f,
                   Y=Y, elbo.vec=elbo.vec)
  return(out.list)
}


flashImpute <- function(Y, idx.backfit=TRUE, greedy.Kmax=50, backfit.maxiter=500){

  out.list <- greedy.init(Y, greedy.Kmax=greedy.Kmax)
  if(idx.backfit==TRUE){out.list <- backfit(Y, out.list, maxiter=backfit.maxiter)}
  return(out.list)
}



# compile functions for a modest speed boost
################################################################################################
library(compiler)

backfit <- cmpfun(backfit)
ELBO <- cmpfun(ELBO)
flashImpute <- cmpfun(flashImpute)
greedy.init <- cmpfun(greedy.init)
init.K <- cmpfun(init.K)
SOL.TAU <- cmpfun(SOL.TAU)
SOLQL.R <- cmpfun(SOLQL.R)

