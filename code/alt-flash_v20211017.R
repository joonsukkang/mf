
library(ebnm)

# run all n updates at the same time (sol.reg)


# main function
alt.flash <- function(fit.init,
                      maxiter=1000,
                      g.common=FALSE,
                      print.elbo=FALSE,
                      sol.reg.max.iter=5,
                      order.output=FALSE
){

  time.vec <- Sys.time()
  elbo.vec <- c()

  for(i in 1:length(fit.init)) assign(names(fit.init)[i], fit.init[[i]])
  g.l.pf <- G.PF(g.l)
  g.f.pf <- G.PF(g.f)

  elbo.tol <- sqrt(.Machine$double.eps) *prod(dim(Y))

  # update g
  if(g.common==TRUE){
    g.lf <- SOL.G.COMMON(g.l.pf, Y, A.l, A.f, B.l, B.f, tau)
    g.l <- g.lf; g.f <- g.lf
  }
  if(g.common==FALSE){
    g.l <- SOL.G(g.l.pf,    Y, A.l, A.f, B.f, tau)
    g.f <- SOL.G(g.f.pf, t(Y), A.f, A.l, B.l, tau)
  }

  # updates
  for (i.iter in 1:maxiter){

    # update q
    obj.temp.l <- SOL.Q(Y, A.l, A.f, B.f, tau, g.l, g.l.pf, sol.reg.max.iter=sol.reg.max.iter)
    A.l <- obj.temp.l$A.l;  B.l <- obj.temp.l$B.l;  C.l <- obj.temp.l$C.l;  KL.l <- obj.temp.l$KL

    # check NULL factor
    k.null <- diag(B.l)==0
    if(sum(k.null)>0){
      print(paste0('drop ', sum(k.null), ' factor(s); ', sum(!k.null), ' factor(s) remaining'))

      A.l <- A.l[,!k.null, drop=FALSE]; C.l <- C.l[,!k.null, drop=FALSE]; B.l <- B.l[!k.null, !k.null, drop=FALSE]
      A.f <- A.f[,!k.null, drop=FALSE]; C.f <- C.f[,!k.null, drop=FALSE]; B.f <- B.f[!k.null, !k.null, drop=FALSE]

      g.l <- g.l[!k.null]; g.f <- g.f[!k.null]

      # refit q
      obj.temp.l <- SOL.Q(Y, A.l, A.f, B.f, tau, g.l, g.l.pf, sol.reg.max.iter=sol.reg.max.iter)
      A.l <- obj.temp.l$A.l;  B.l <- obj.temp.l$B.l;  C.l <- obj.temp.l$C.l;  KL.l <- obj.temp.l$KL

      obj.temp.f <- SOL.Q(t(Y), A.f, A.l, B.l, tau, g.f, g.f.pf, sol.reg.max.iter=sol.reg.max.iter)
      A.f <- obj.temp.f$A.l;  B.f <- obj.temp.f$B.l;  C.f <- obj.temp.f$C.l;  KL.f <- obj.temp.f$KL
    }


    obj.temp.f <- SOL.Q(t(Y), A.f, A.l, B.l, tau, g.f, g.f.pf, sol.reg.max.iter=sol.reg.max.iter)
    A.f <- obj.temp.f$A.l;  B.f <- obj.temp.f$B.l;  C.f <- obj.temp.f$C.l;  KL.f <- obj.temp.f$KL


    # check NULL factor
    k.null <- diag(B.f)==0
    if(sum(k.null)>0){
      print(paste0('drop ', sum(k.null), ' factor(s); ', sum(!k.null), ' factor(s) remaining'))

      A.l <- A.l[,!k.null, drop=FALSE]; C.l <- C.l[,!k.null, drop=FALSE]; B.l <- B.l[!k.null, !k.null, drop=FALSE]
      A.f <- A.f[,!k.null, drop=FALSE]; C.f <- C.f[,!k.null, drop=FALSE]; B.f <- B.f[!k.null, !k.null, drop=FALSE]

      g.l <- g.l[!k.null]; g.f <- g.f[!k.null]

      # refit q
      obj.temp.l <- SOL.Q(Y, A.l, A.f, B.f, tau, g.l, g.l.pf, sol.reg.max.iter=sol.reg.max.iter)
      A.l <- obj.temp.l$A.l;  B.l <- obj.temp.l$B.l;  C.l <- obj.temp.l$C.l;  KL.l <- obj.temp.l$KL

      obj.temp.f <- SOL.Q(t(Y), A.f, A.l, B.l, tau, g.f, g.f.pf, sol.reg.max.iter=sol.reg.max.iter)
      A.f <- obj.temp.f$A.l;  B.f <- obj.temp.f$B.l;  C.f <- obj.temp.f$C.l;  KL.f <- obj.temp.f$KL
    }


    # compute elbo
    elbo <- ELBO(Y, trYTY, tau, A.l, A.f, B.l, B.f, KL.l, KL.f)
    if(print.elbo==TRUE){print(paste0(Sys.time(), "; iter=", i.iter, "; elbo=", elbo))}

    # record time and elbo
    time.vec <- c(time.vec, Sys.time())
    elbo.vec <- c(elbo.vec, elbo)

    # end loop if converged
    if(i.iter>1){ if( elbo-elbo.old < elbo.tol) break }
    elbo.old <- elbo

    # update tau
    tau <- SOL.TAU(trYTY=trYTY, Y, A.l, A.f, B.l, B.f)

    # update g
    if(g.common==TRUE){
      g.lf <- SOL.G.COMMON(g.l.pf, Y, A.l, A.f, B.l, B.f, tau)
      g.l <- g.lf; g.f <- g.lf
    }
    if(g.common==FALSE){
      g.l <- SOL.G(g.l.pf,    Y, A.l, A.f, B.f, tau)
      g.f <- SOL.G(g.f.pf, t(Y), A.f, A.l, B.l, tau)
    }

  }

  if(i.iter==maxiter){ print("Maximum number of iterations reached.") }
  print(paste0('backfit completed; elbo=', elbo))

  # scale and order factors
  scales <- sqrt(diag(B.l) * diag(B.f))
  if(order.output==TRUE){
    ordering <- order(scales, decreasing=TRUE)
    scales <- scales[ordering]
    A.l <- A.l[,ordering,drop=FALSE]; C.l <- C.l[,ordering,drop=FALSE]
    A.f <- A.f[,ordering,drop=FALSE]; C.f <- C.f[,ordering,drop=FALSE]

    B.l <- crossprod(A.l); diag(B.l) <- colSums(C.l)
    B.f <- crossprod(A.f); diag(B.f) <- colSums(C.f)

    g.l <- g.l[ordering]; g.f <- g.f[ordering]
  }

  time.vec <- as.numeric(time.vec[-1] - time.vec[1])

  out.list <- list(A.l=A.l, A.f=A.f,
                   B.l=B.l, B.f=B.f,
                   C.l=C.l, C.f=C.f,
                   g.l=g.l, g.f=g.f,
                   tau=tau,
                   elbo.vec=elbo.vec,
                   time.vec=time.vec,
                   elbo.tol=elbo.tol,
                   scales=scales)

  return(out.list)
}

# return nonnegative direction estiamtes
INIT.L.NonNeg <- function(Y, verbose=1){

  max.K <- 20
  tol.scale.ratio <- 1e-4

  n <- nrow(Y)

  # residual Y
  Yres <- Y

  # root factor
  u <- eigen(Yres)$vectors[,1]
  u <- u * sign(u[which.max(abs(u))])
  u[u<0] <- 0
  u <- u/sqrt(sum(u^2))
  L <- matrix(u, ncol=1)
  d <- scaleL.qp(Y, L, scale.min=0, print.scale=FALSE)$scale
  Yres <- Y - L %*% diag(d, nrow=length(d)) %*% t(L)

  # add factors
  for (k in 1:max.K){
    ## find direction
    u <- eigen(Yres)$vectors[,1]
    u <- u * sign(u[which.max(abs(u))])
    u[u<0] <- 0
    u <- u/sqrt(sum(u^2))

    ## find scale
    L <- cbind(L, matrix(u, ncol=1))
    d <- scaleL.qp(Y, L, scale.min=0, print.scale=FALSE)$scale
    if(verbose==1){print(paste0('min/max sacle ratio: ', min(d)/max(d)))}
    if(min(d)/max(d) < tol.scale.ratio) {
      L <- L[,-k]
      d <- scaleL.qp(Y, L, scale.min=0, print.scale=FALSE)$scale
      break
    }
    Yres <- Y - L %*% diag(d, nrow=length(d)) %*% t(L)
  }

  if(verbose==1){print(paste0(ncol(L), ' factors fitted'))}
  out.list <- list(L=L)
  return(out.list)
}


# given 'init.L' and 'Y', create 'init.fit' object
INIT.L <- function(init.L, Y){
  n <- nrow(Y); p <- ncol(Y); K <- ncol(init.L)
  A.l <- init.L; A.f <- init.L

  g.lf <- list()
  for (k in 1:K){g.lf[[k]] <- gammamix(pi=c(mean(A.l[,k]==0), mean(A.l[,k]!=0)),
                                       shape=c(1,1), scale=c(0, mean(A.l[A.l[,k]!=0,k])))}
  g.l <- g.lf; g.f <- g.lf

  trYTY <- sum(diag(Y %*% t(Y)))
  C.l <- A.l^2; B.l <- crossprod(A.l)
  C.f <- A.f^2; B.f <- crossprod(A.f)
  tau <- SOL.TAU(trYTY=trYTY, Y, A.l, A.f, B.l, B.f)

  out.list <- list(Y=Y, trYTY=trYTY,
                   n=n, p=p, K=K,
                   A.l=A.l, A.f=A.f,
                   B.l=B.l, B.f=B.f,
                   C.l=C.l, C.f=C.f,
                   tau=tau,
                   g.l=g.l, g.f=g.f
  )
  return(out.list)
}




# initialization from fit.flashier; scale applied to F
INIT.flashier <- function(fit.init){

  Y <- fit.init$flash.fit$Y
  K <- fit.init$n.factors
  trYTY <- fit.init$flash.fit$Y2 # = sum(diag(crossprod(Y,Y)))
  n <- nrow(Y)
  p <- ncol(Y)

  A.l <- fit.init$flash.fit$EF[[1]]
  A.f <- fit.init$flash.fit$EF[[2]]
  C.l <- fit.init$flash.fit$EF2[[1]]
  C.f <- fit.init$flash.fit$EF2[[2]]

  B.l <- crossprod(A.l, A.l)
  B.f <- crossprod(A.f, A.f)
  diag(B.l) <- colSums(C.l)
  diag(B.f) <- colSums(C.f)

  tau <- fit.init$flash.fit$tau
  g.l <- fit.init$fitted.g[[1]]
  g.f <- fit.init$fitted.g[[2]]

  out.list <- list(Y=Y, trYTY=trYTY,
                   n=n, p=p, K=K,
                   A.l=A.l, A.f=A.f,
                   B.l=B.l, B.f=B.f,
                   C.l=C.l, C.f=C.f,
                   tau=tau,
                   g.l=g.l, g.f=g.f
  )
  return(out.list)
}


# return the "prior_family" (for the ebnm argument)
G.PF <- function(g){

  if( (prod(sapply(g, function(x) class(x)=="normalmix"))==1) ){
    if (prod(sapply(g, function(x) length(x$sd)==1))==1) {g.pf <- 'normal'} }

  if( (prod(sapply(g, function(x) class(x)=="unimix"))==1) ){
    if ( prod(sapply(sapply(g, function(x) x$a>=0), prod))==1 &
         prod(sapply(sapply(g, function(x) x$b>=0), prod))==1  ) {g.pf <- 'unimodal_nonnegative'} }

  if( (prod(sapply(g, function(x) class(x)=="gammamix"))==1) ){
    if ( prod(sapply(g, function(x) x$shape == c(1,1)))==1 &
         prod(sapply(g, function(x) x$scale[1] == 0))==1) {g.pf <- 'point_exponential'} }

  return(g.pf)
}

####################### WORKING #######################


# coordinate descent algorithm for updating q_l
SOL.Q <- function(Y, A.l, A.f, B.f, tau, g.l, g.l.pf, sol.reg.max.iter){

  n <- nrow(Y)
  K <- length(g.l)
  KL <- rep(0, K)
  C.l <- matrix(0, nrow=n, ncol=K)

  for (it in 1:sol.reg.max.iter){
    for (k in 1:K){

      x <- (Y %*% A.f[,k] - A.l[,-k] %*% B.f[-k,k])/B.f[k,k]
      s <- rep(1/sqrt(tau*B.f[k,k]), length(x))

      temp.obj <- ebnm(x=x, s=s, prior_family=g.l.pf, fix_g=TRUE, g_init=g.l[[k]],
                       output=c('posterior_mean', "posterior_second_moment", "log_likelihood"))

      Et <- A.l[,k] <- temp.obj$posterior$mean
      Et2 <- C.l[,k] <- temp.obj$posterior$second_moment
      KL[k] <- temp.obj$log_likelihood + sum(log(2*pi*s^2)/2) + sum((x^2 - 2*x*Et + Et2)/(2*s^2))
    }
  }

  B.l <- crossprod(A.l)
  diag(B.l) <- colSums(C.l)

  out.list <- list(A.l=A.l, B.l=B.l, C.l=C.l, KL=sum(KL))
  return(out.list)
}






# update g_l (or g_f)
SOL.G <- function(g.l.pf, Y, A.l, A.f, B.f, tau){

  K <- ncol(A.l)
  g.l <- list()

  for (k in 1:K){

    x <- c(Y%*%A.f[,k] - A.l[,-k, drop=FALSE] %*% B.f[k,-k])/B.f[k,k]
    s <- rep(1/sqrt(tau*B.f[k,k]), length(x))

    obj.temp <- ebnm(x=x, s=s, prior_family=g.l.pf, fix_g=FALSE)
    g.l[[k]]<- obj.temp$fitted_g
  }

  return(g.l)
}


SOL.G.COMMON <- function(g.l.pf, Y, A.l, A.f, B.l, B.f, tau){

  K <- ncol(A.l)
  g.lf <- list()

  for (k in 1:K){
    x.l <- c(Y%*%A.f[,k] - A.l[,-k, drop=FALSE] %*% B.f[k,-k])/B.f[k,k]
    s.l <- rep(1/sqrt(tau*B.f[k,k]), length(x.l))
    x.f <- c(t(Y)%*%A.l[,k] - A.f[,-k, drop=FALSE] %*% B.l[k,-k])/B.l[k,k]
    s.f <- rep(1/sqrt(tau*B.l[k,k]), length(x.f))

    obj.temp <- ebnm(x=c(x.l, x.f), s=c(s.l, s.f), prior_family=g.l.pf, fix_g=FALSE)
    g.lf[[k]] <- obj.temp$fitted_g
  }

  return(g.lf)
}






# update tau
SOL.TAU <- function(trYTY, Y, A.l, A.f, B.l, B.f){

  n <- nrow(Y)
  p <- ncol(Y)

  if(n<p){  tau <- ((trYTY - 2 * sum(diag(tcrossprod(crossprod(Y, A.l),A.f))) + sum(B.l*B.f))/(n*p))^(-1)}
  if(n>=p){ tau <- ((trYTY - 2 * sum(diag(tcrossprod(tcrossprod(t(A.f), Y), t(A.l)))) + sum(B.l*B.f))/(n*p))^(-1)}
  return(tau)
}


# compute elbo
ELBO <- function(Y, trYTY, tau,
                 A.l, A.f, B.l, B.f,
                 KL.l, KL.f){
  n <- nrow(Y)
  p <- ncol(Y)
  if(n<p){  chunk <- (trYTY - 2 * sum(diag(tcrossprod(crossprod(Y, A.l),A.f))) + sum(B.l*B.f))}
  if(n>=p){ chunk <- (trYTY - 2 * sum(diag(tcrossprod(tcrossprod(t(A.f), Y), t(A.l)))) + sum(B.l*B.f))}
  elbo.0 <- -(n*p/2)*log(2*pi) + (n*p/2)*log(tau) - (tau/2)*chunk

  elbo <- elbo.0 + KL.l + KL.f
  return(elbo)
}





# compile functions for a modest speed boost
################################################################################################
library(compiler)

alt.flash <- cmpfun(alt.flash)
ELBO <- cmpfun(ELBO)
G.PF <- cmpfun(G.PF)
INIT.flashier <- cmpfun(INIT.flashier)
INIT.L <- cmpfun(INIT.L)
SOL.G <- cmpfun(SOL.G)
SOL.Q <- cmpfun(SOL.Q)
SOL.TAU <- cmpfun(SOL.TAU)

