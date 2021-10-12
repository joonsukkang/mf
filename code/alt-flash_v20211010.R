
library(ebnm)

# main function
alt.flash <- function(flashier.fit=NULL, # initialize with flashier obj
                      Y=NULL, # initialize with (Y,K,prior.family)
                      g.l=NULL,
                      g.f=NULL,
                      seed=NULL, # seed for random initialization
                      maxiter=1000,
                      g.common=FALSE
                      ){

  time.vec <- Sys.time()
  elbo.vec <- c()

  # initialize from flashier fit
  if(!is.null(flashier.fit)){ obj.temp <- INIT.flashier(flashier.fit) }
  # random initialization
  if(is.null(flashier.fit)){ obj.temp <- INIT.rand(Y, g.l, g.f, seed) }

  for(i in 1:length(obj.temp)) assign(names(obj.temp)[i], obj.temp[[i]])
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
    obj.temp.l <- SOL.Q(Y, A.l, A.f, B.f, tau, g.l, g.l.pf)
    A.l <- obj.temp.l$A.l;  B.l <- obj.temp.l$B.l;  C.l <- obj.temp.l$C.l;  KL.l <- obj.temp.l$KL

        # check NULL factor
        k.null <- diag(B.l)==0
        if(sum(k.null)>0){
          print(paste0('drop ', sum(k.null), ' factor(s); ', sum(!k.null), ' factor(s) remaining'))

          A.l <- A.l[,!k.null, drop=FALSE]; C.l <- C.l[,!k.null, drop=FALSE]; B.l <- B.l[!k.null, !k.null, drop=FALSE]
          A.f <- A.f[,!k.null, drop=FALSE]; C.f <- C.f[,!k.null, drop=FALSE]; B.f <- B.f[!k.null, !k.null, drop=FALSE]

          g.l <- g.l[!k.null]; g.f <- g.f[!k.null]

          # refit q
          obj.temp.l <- SOL.Q(Y, A.l, A.f, B.f, tau, g.l, g.l.pf)
          A.l <- obj.temp.l$A.l;  B.l <- obj.temp.l$B.l;  C.l <- obj.temp.l$C.l;  KL.l <- obj.temp.l$KL

          obj.temp.f <- SOL.Q(t(Y), A.f, A.l, B.l, tau, g.f, g.f.pf)
          A.f <- obj.temp.f$A.l;  B.f <- obj.temp.f$B.l;  C.f <- obj.temp.f$C.l;  KL.f <- obj.temp.f$KL
        }


    obj.temp.f <- SOL.Q(t(Y), A.f, A.l, B.l, tau, g.f, g.f.pf)
    A.f <- obj.temp.f$A.l;  B.f <- obj.temp.f$B.l;  C.f <- obj.temp.f$C.l;  KL.f <- obj.temp.f$KL


        # check NULL factor
        k.null <- diag(B.f)==0
        if(sum(k.null)>0){
          print(paste0('drop ', sum(k.null), ' factor(s); ', sum(!k.null), ' factor(s) remaining'))

          A.l <- A.l[,!k.null, drop=FALSE]; C.l <- C.l[,!k.null, drop=FALSE]; B.l <- B.l[!k.null, !k.null, drop=FALSE]
          A.f <- A.f[,!k.null, drop=FALSE]; C.f <- C.f[,!k.null, drop=FALSE]; B.f <- B.f[!k.null, !k.null, drop=FALSE]

          g.l <- g.l[!k.null]; g.f <- g.f[!k.null]

          # refit q
          obj.temp.l <- SOL.Q(Y, A.l, A.f, B.f, tau, g.l, g.l.pf)
          A.l <- obj.temp.l$A.l;  B.l <- obj.temp.l$B.l;  C.l <- obj.temp.l$C.l;  KL.l <- obj.temp.l$KL

          obj.temp.f <- SOL.Q(t(Y), A.f, A.l, B.l, tau, g.f, g.f.pf)
          A.f <- obj.temp.f$A.l;  B.f <- obj.temp.f$B.l;  C.f <- obj.temp.f$C.l;  KL.f <- obj.temp.f$KL
        }


    # compute elbo
    elbo <- ELBO(Y, trYTY, tau, A.l, A.f, B.l, B.f, KL.l, KL.f)
    #print(elbo)

    # record time and elbo
    time.vec <- c(time.vec, Sys.time())
    elbo.vec <- c(elbo.vec, elbo)

    # end loop if converged
    if(i.iter>1){ if(abs(elbo-elbo.old)< elbo.tol) break }
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
  ordering <- order(scales, decreasing=TRUE)
  scales <- scales[ordering]
  A.l <- A.l[,ordering,drop=FALSE]; C.l <- C.l[,ordering,drop=FALSE]
  A.f <- A.f[,ordering,drop=FALSE]; C.f <- C.f[,ordering,drop=FALSE]

  B.l <- crossprod(A.l); diag(B.l) <- colSums(C.l)
  B.f <- crossprod(A.f); diag(B.f) <- colSums(C.f)

  g.l <- g.l[ordering]; g.f <- g.f[ordering]


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


# sample from point exponential (pe)
sample.pe <- function(g, nsamp){

  z <- rbinom(n=nsamp, size=1, prob=g$pi[2])
  x <- rgamma(n=nsamp, shape=g$shape[2], scale=g$scale[2])
  x[z==0] <- 0
  return(x)
}




# random initialization
INIT.rand <- function(Y, g.l, g.f, elbo.tol=1, max.update=30, min.update=10, seed=NULL){

  if(!is.null(seed)){set.seed(seed)}

  g.l.pf <- G.PF(g.l)
  g.f.pf <- G.PF(g.f)
  if(g.l.pf=='point_exponential'){sample.L <- sample.pe}
  if(g.f.pf=='point_exponential'){sample.F <- sample.pe}

  n <- nrow(Y); p <- ncol(Y); K <- length(g.l)
  trYTY <- sum(diag(crossprod(Y,Y)))

  A.l <- matrix(0, nrow=n, ncol=K)
  A.f <- matrix(0, nrow=p, ncol=K)
  for (k in 1:K){
    A.l[,k] <- sample.L(g.l[[k]], nsamp=n)
    A.f[,k] <- sample.F(g.f[[k]], nsamp=p)
  }
  B.f <- crossprod(A.f); B.l <- crossprod(A.l)
  C.f <- A.f^2; C.l <- A.l^2
  tau <- 1/sd(Y)

  # update with g.l and g.f fixed
  # until elbo stabilizes
  for (i in 1:max.update){
    obj.temp.l <- SOL.Q(Y, A.l, A.f, B.f, tau, g.l, g.l.pf)
    A.l <- obj.temp.l$A.l;  B.l <- obj.temp.l$B.l;  C.l <- obj.temp.l$C.l;  KL.l <- obj.temp.l$KL

    obj.temp.f <- SOL.Q(t(Y), A.f, A.l, B.l, tau, g.f, g.f.pf)
    A.f <- obj.temp.f$A.l;  B.f <- obj.temp.f$B.l;  C.f <- obj.temp.f$C.l;  KL.f <- obj.temp.f$KL

    tau <- SOL.TAU(trYTY=trYTY, Y, A.l, A.f, B.l, B.f)
    elbo <- ELBO(Y, trYTY, tau, A.l, A.f, B.l, B.f, KL.l, KL.f)
    #print(elbo)

    if(i > min.update){if(abs(elbo.old-elbo) < elbo.tol) break}
    elbo.old <- elbo
  }
  print(paste0('initialization completed; elbo=', elbo))

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

# returns a coordinate-wise solution
SOL.C <- function(y, A.f.k, B.f.kk, B.f.minus.k, b.minus.k, tau, g, g.pf){

  x <- (sum(y*A.f.k)-sum(B.f.minus.k*b.minus.k))/B.f.kk
  s <- rep(1/sqrt(tau*B.f.kk), length(x))

  temp.obj <- ebnm(x=x, s=s, prior_family=g.pf, fix_g=TRUE, g_init=g)

  Et <- temp.obj$posterior$mean
  Et2 <- temp.obj$posterior$sd^2 + Et^2
  KL <- temp.obj$log_likelihood + sum(log(2*pi*s^2)/2) + sum((x^2 - 2*x*Et + Et2)/(2*s^2))

  out.list <- list(Et=Et, Et2=Et2, KL=KL)
  return(out.list)
}



# returns a coordinate descent solution for a regression problem
SOL.REG <- function(y, b, A.f, B.f, tau, g.l, g.l.pf, max.iter=100, tol=1e-10){

  K <- length(g.l)
  b.old <- b
  for (i.iter in 1:max.iter){
    for (k in 1:K){
      b[k] <- SOL.C(y=y,
                    A.f.k=A.f[,k],
                    B.f.kk=B.f[k,k],
                    B.f.minus.k=B.f[k,-k],
                    b.minus.k = b[-k],
                    tau = tau,
                    g = g.l[[k]],
                    g.pf = g.l.pf)$Et
    }
    if(sum((b-b.old)^2) < tol) break
    b.old <- b
  }

  Et <- rep(0, K)
  Et2 <- rep(0,K)
  KL <- rep(0,K)
  for (k in 1:K){
    temp.obj <- SOL.C(y=y,
                  A.f.k=A.f[,k],
                  B.f.kk=B.f[k,k],
                  B.f.minus.k=B.f[k,-k],
                  b.minus.k = b[-k],
                  tau = tau,
                  g = g.l[[k]],
                  g.pf = g.l.pf)
    Et[k]  <- temp.obj$Et
    Et2[k] <- temp.obj$Et2
    KL[k]  <- temp.obj$KL
  }

  out.list <- list(Et=Et, Et2=Et2, KL=KL)
  return(out.list)
}


# coordinate descent algorithm for updating q_l
SOL.Q <- function(Y, A.l, A.f, B.f, tau, g.l, g.l.pf){

  n <- nrow(Y)
  K <- length(g.l)
  KL <- 0
  C.l <- matrix(0, nrow=n, ncol=K)

  for (i in 1:n){
    temp.obj <- SOL.REG(y=Y[i,], b=A.l[i,], A.f, B.f, tau, g.l, g.l.pf)
    A.l[i,] <- temp.obj$Et
    C.l[i,] <- temp.obj$Et2
    KL <- KL + sum(temp.obj$KL)
  }

  B.l <- crossprod(A.l)
  diag(B.l) <- colSums(C.l)

  out.list <- list(A.l=A.l, B.l=B.l, C.l=C.l, KL=KL)
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
INIT.rand <- cmpfun(INIT.rand)
SOL.C <- cmpfun(SOL.C)
SOL.G <- cmpfun(SOL.G)
SOL.Q <- cmpfun(SOL.Q)
SOL.REG <- cmpfun(SOL.REG)
SOL.TAU <- cmpfun(SOL.TAU)

