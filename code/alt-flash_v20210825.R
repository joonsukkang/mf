

# main function
alt.flash <- function(fit.init,
                      n.cores=1,
                      n.jobs=NULL,
                      maxiter.mean.sol.c=0,
                      maxiter=1000,
                      drop.complete=TRUE,
                      when.normal.closed.form=TRUE){

  if(is.null(n.jobs)==TRUE){n.jobs <- n.cores}
  time.vec <- Sys.time()
  elbo.vec <- c()

  # initialize from flashier fit
  obj.temp <- INIT(fit.init)
  for(i in 1:length(obj.temp)) assign(names(obj.temp)[i], obj.temp[[i]])
  g.l.pf <- G.PF(g.l)
  g.f.pf <- G.PF(g.f)

  elbo.tol <- sqrt(.Machine$double.eps) *prod(dim(Y))


  # updates
  for (i.iter in 1:maxiter){

    # update q
    obj.temp.l <- SOL.Q(  Y,  A.l, A.f, B.f, tau, g.l, g.l.pf, maxiter.mean.sol.c, n.cores, n.jobs, drop.complete, when.normal.closed.form)
    A.l <- obj.temp.l$A;  B.l <- obj.temp.l$B;  C.l <- obj.temp.l$C;  KL.l <- obj.temp.l$KL

    obj.temp.f <- SOL.Q(t(Y), A.f, A.l, B.l, tau, g.f, g.f.pf, maxiter.mean.sol.c, n.cores, n.jobs, drop.complete, when.normal.closed.form)
    A.f <- obj.temp.f$A;  B.f <- obj.temp.f$B;  C.f <- obj.temp.f$C;  KL.f <- obj.temp.f$KL

    # compute elbo
    elbo <- ELBO(Y, trYTY, tau, A.l, A.f, B.l, B.f, KL.l, KL.f)

    # record time and elbo
    time.vec <- c(time.vec, Sys.time())
    elbo.vec <- c(elbo.vec, elbo)

    # end loop if converged
    if(i.iter>1){ if(abs(elbo-elbo.old)< elbo.tol) break }
    elbo.old <- elbo

    # update tau
    tau <- SOL.TAU(trYTY=trYTY, Y, A.l, A.f, B.l, B.f)

    # update g
    g.l <- SOL.G(g.l, g.l.pf,   Y,  A.l, A.f, B.f, C.l, tau, obj.temp.l$mean.w)
    g.f <- SOL.G(g.f, g.f.pf, t(Y), A.f, A.l, B.l, C.f, tau, obj.temp.f$mean.w)
  }

  if(i.iter==maxiter){ print("Maximum number of iterations reached.") }

  time.vec <- as.numeric(time.vec[-1] - time.vec[1])

  out.list <- list(A.l=A.l, A.f=A.f,
                   B.l=B.l, B.f=B.f,
                   C.l=C.l, C.f=C.f,
                   g.l=g.l, g.f=g.f,
                   tau=tau,
                   elbo.vec=elbo.vec,
                   time.vec=time.vec,
                   elbo.tol=elbo.tol)

  return(out.list)
}



# initialization from fit.flashier; scale applied to F
INIT <- function(fit.init){

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
  return(g.pf)
}

####################### WORKING #######################


# solve a regression problem for a batch
SOL.REG.B <- function(mat.y, mat.b, B.f, tau, g.l, g.l.pf, maxiter.mean.sol.c, drop.complete){

  bs <- nrow(mat.b) # block size
  K <- ncol(mat.b)

  i.iter <- 0
  idx.included <- 1:bs

  while( (i.iter < maxiter.mean.sol.c ) & (length(idx.included) > 0) ){

    i.iter <- i.iter + 1
    temp.mat.b <- mat.b[idx.included, ,drop=FALSE]
    old.temp.mat.b <- temp.mat.b

    for (k in 1:K){
      temp.mat.b[, k] <- MEAN.SOL.B(chunk = c(mat.y[idx.included,k]-tcrossprod(temp.mat.b[,-k], B.f[k,-k, drop=FALSE])),
                                    Bkk=B.f[k,k],
                                    tau=tau,
                                    g=g.l[[k]],
                                    g.pf=g.l.pf)
    }
    mat.b[idx.included,] <- temp.mat.b # save results
    if(drop.complete){ # drop completed row (if drop.complete option is selected)
      idx.included <- idx.included[rowSums((old.temp.mat.b - temp.mat.b)^2)>1e-6]
    }
  }

  mat.Et2 <- matrix(0, nrow=bs, ncol=K)
  KL <- 0
  list.w <- list()

  for (k in 1:K){
    obj.temp <- SOL.B(chunk=as.vector(mat.y[,k]- tcrossprod(mat.b[,-k],  B.f[k,-k, drop=FALSE])),
                      Bkk=B.f[k,k],
                      tau=tau, g=g.l[[k]], g.pf=g.l.pf)
    mat.b[,k] <- obj.temp$Et
    mat.Et2[,k] <- obj.temp$Et2
    KL <- KL + obj.temp$KL
    if(g.l.pf=='unimodal_nonnegative'){ list.w[[k]] <- rowSums(obj.temp$w)}
  }

  out.list <- list(Et=mat.b, Et2=mat.Et2, KL=KL)
  if(g.l.pf=='unimodal_nonnegative'){out.list <- list(Et=mat.b, Et2=mat.Et2, KL=KL, list.w=list.w)}

  return(out.list)
}



# return posterior mean for one batch
MEAN.SOL.B <- function(chunk, Bkk, tau, g, g.pf){

  x <- chunk/Bkk
  s <- rep(1/sqrt(tau*Bkk), length(x))
  data <- set_data(betahat=x, sebetahat=s)
  Et <- postmean(g, data)

  return(Et)
}

# return posterior mean, mean2, and KL for one batch
SOL.B <- function(chunk, Bkk, tau, g, g.pf){

  x <- chunk/Bkk
  s <- rep(1/sqrt(tau*Bkk), length(x))
  data <- set_data(betahat=x, sebetahat=s)

  if(g.pf=='normal'){
    Et <- postmean(g, data)
    Et2 <- postmean2(g, data)
    KL <-  calc_loglik(g, data) + sum(log(2*pi*s^2)/2) + sum((x^2 - 2*x*Et + Et2)/(2*s^2))
    out.list <- list(Et=Et, Et2=Et2, KL=KL)
  }
  if(g.pf=='unimodal_nonnegative'){
    w <- comp_postprob(g, data)
    Et <- colSums(w * comp_postmean(g, data))
    Et2 <- colSums(w * comp_postmean2(g, data)); Et2[Et2<0] <- 0
    KL <-  calc_loglik(g, data) + sum(log(2*pi*s^2)/2) + sum((x^2 - 2*x*Et + Et2)/(2*s^2))
    out.list <- list(Et=Et, Et2=Et2, KL=KL, w=w)
  }

  return(out.list)
}


# update q_l (or q_f)
SOL.Q <- function(Y, A.l, A.f, B.f, tau, g.l, g.l.pf, maxiter.mean.sol.c, n.cores, n.jobs, drop.complete, when.normal.closed.form){

  n <- nrow(Y)
  K <- ncol(A.l)
  C.l <- matrix(0, nrow=n, ncol=K)
  KL <- 0

  if( (g.l.pf=='normal' & when.normal.closed.form==FALSE) | g.l.pf=='unimodal_nonnegative'){

    Y.A.f <- tcrossprod(Y, t(A.f))

    # split into 'n.jobs' jobs
    bite.size <- floor(n/n.jobs)
    SOL.REG.i <- function(i){
      idx.from <- (i-1)*bite.size + 1
      idx.to <- i*bite.size; if(i==n.jobs){idx.to <- n}
      obj.temp <- SOL.REG.B(mat.y=Y.A.f[idx.from:idx.to, , drop=FALSE], mat.b=A.l[idx.from:idx.to, , drop=FALSE],
                            B.f=B.f, tau=tau, g.l=g.l, g.l.pf=g.l.pf,
                            maxiter.mean.sol.c=maxiter.mean.sol.c, drop.complete=drop.complete)
      return(obj.temp)
    }

    results <- mclapply(1:n.jobs, SOL.REG.i, mc.cores=n.cores)

    A.l <- do.call(rbind, lapply(results, function(x) x$Et))
    C.l <- do.call(rbind, lapply(results, function(x) x$Et2))
    KL <- sum(unlist(lapply(results, function(x) x$KL)))

    if(g.l.pf=='unimodal_nonnegative'){
      mean.w <- list()
      for (k in 1:K){
        mean.w[[k]] <- colSums(do.call(rbind, lapply(results, function(x,k) x$list.w[[k]], k)))
        mean.w[[k]] <- mean.w[[k]]/sum(mean.w[[k]])
      }
    }

  }

  if(g.l.pf=='normal' & when.normal.closed.form==TRUE){ # use closed-form solution
    tau.l <- sapply(g.l, function(x) x$sd)^(-2)
    A.l <- tcrossprod(tcrossprod(Y, t(A.f)),  t(solve(B.f+diag(tau.l/tau))))

    for (k in 1:K){
      tau.k <- g.l[[k]]$sd^(-2)
      Et <- A.l[,k]
      Et2 <- Et^2 + 1/(tau.k + tau*B.f[k,k])
      KL <- KL -(1/2) * ( tau.k * sum(Et2) - n -n*log(tau.k) + n*log(tau.k+tau*B.f[k,k]))
      C.l[,k] <- Et2
    }
  }

  B.l <- crossprod(A.l, A.l); diag(B.l) <- colSums(C.l)

  if(g.l.pf=='normal'){out.list = list(A=A.l, B=B.l, C=C.l, KL=KL) }
  if(g.l.pf=='unimodal_nonnegative'){ out.list = list(A=A.l, B=B.l, C=C.l, KL=KL, mean.w=mean.w) }
  return(out.list)
}


# update g_l (or g_f)
SOL.G <- function(g.l, g.l.pf, Y, A.l, A.f, B.f, C.l, tau, mean.w=NULL){

  n <- nrow(Y)
  K <- ncol(A.l)

  if(g.l.pf=='normal'){
    for (k in 1:K){ g.l[[k]]$sd <- mean(C.l[,k])^(1/2) }
    }

  if(g.l.pf=='unimodal_nonnegative'){
    for (k in 1:K){ g.l[[k]]$pi <-mean.w[[k]] }
  }

  return(g.l)
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




################################################################################################
# modified version of flashier;:flash.backfit

flash.backfit.e <- function (flash, kset = NULL, method = c("extrapolate", "sequential",
                                                                        "dropout", "random", "parallel"), warmstart = TRUE, conv.crit.fn = calc.obj.diff,
                                         tol = set.default.tol(flash), maxiter = 500, verbose.lvl = get.verbose.lvl(flash)){

  ################################################ EDITED FOR ELBO-TIME TRACKING
  time.vec <- Sys.time()
  elbo.vec <- c()
  ################################################

  flash <- get.fit(flash)
  if (is.null(kset)) {
    if (get.n.factors(flash) > 0) {
      kset <- 1:get.n.factors(flash)
    }
    else {
      announce.no.backfit(verbose.lvl)
      verbose.lvl <- 0
    }
  }
  else {
    must.be.valid.kset(flash, kset)
  }
  must.be.integer(maxiter, lower = 1, allow.null = FALSE)
  must.be.integer(verbose.lvl, lower = -1, upper = 3, allow.null = FALSE)
  method <- match.arg(method)
  if (method == "parallel") {
    check.parallel.ok(flash, kset)
    if (missing(conv.crit.fn)) {
      conv.crit.fn <- function(new, old, k) {
        return(abs(calc.obj.diff(new, old, k)))
      }
    }
  }
  if (missing(tol)) {
    report.tol.setting(verbose.lvl, tol)
  }
  else {
    must.be.numeric(tol, allow.infinite = FALSE, allow.null = FALSE)
  }
  flash <- set.warmstart(flash, warmstart)
  verbose.fns <- get.verbose.fns(flash)
  verbose.colnames <- get.verbose.colnames(flash)
  verbose.colwidths <- get.verbose.colwidths(flash)
  conv.crit <- rep(Inf, get.n.factors(flash))
  conv.crit[setdiff(1:get.n.factors(flash), kset)] <- 0
  announce.backfit(verbose.lvl, n.factors = length(kset), tol)
  print.table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                     backfit = TRUE)
  if (method == "parallel") {
    kset <- setdiff(kset, which(is.zero(flash)))
    kset <- setdiff(kset, which.k.fixed(flash))
    cl <- parallel::makeCluster(getOption("cl.cores", 2L),
                                type = getOption("cl.type", "PSOCK"), useXDR = FALSE)
  }
  else if (method == "extrapolate") {
    extrapolate.control <- getOption("extrapolate.control",
                                     list())
    extrapolate.param <- set.extrapolate.param(extrapolate.control)
  }
  iter <- 0
  old.obj <- get.obj(flash)
  next.tol.target <- NULL
  if (method == "extrapolate") {
    extrapolate.param <- init.beta(extrapolate.param)
    old.f <- flash
  }
  while (iter < maxiter && max(conv.crit) > tol) {
    iter <- iter + 1
    kset <- get.next.kset(method, kset, conv.crit, tol)
    if (!(method %in% c("parallel", "extrapolate"))) {
      for (k in kset) {
        old.f <- flash
        flash <- update.one.factor(flash, k, iter, verbose.lvl)
        info <- calc.update.info(flash, old.f, conv.crit.fn,
                                 verbose.fns, k)
        conv.crit[k] <- get.conv.crit(info)
        print.table.entry(verbose.lvl, verbose.colwidths,
                          iter, info, k = k, backfit = TRUE)
      }
    }
    else {
      if (method == "parallel") {
        old.f <- flash
        flash <- update.factors.parallel(flash, kset,
                                         cl)
      }
      else if (method == "extrapolate") {
        proposed.f <- extrapolate.f(flash, old.f, extrapolate.param)
        proposed.f <- update.factors.in.kset(proposed.f,
                                             kset)
        old.f <- flash
        if (get.obj(proposed.f) - get.obj(flash) < tol) {
          flash <- update.factors.in.kset(flash, kset)
          extrapolate.param <- decelerate(extrapolate.param)
        }
        else {
          flash <- proposed.f
          extrapolate.param <- accelerate(extrapolate.param)
        }
      }
      info <- calc.update.info(flash, old.f, conv.crit.fn,
                               verbose.fns)
      conv.crit <- get.conv.crit(info)
      print.table.entry(verbose.lvl, verbose.colwidths,
                        iter, info, k = "all", backfit = TRUE)
    }
    if (is.null(next.tol.target) && max(conv.crit) > 0 &&
        max(conv.crit) < Inf) {
      next.tol.target <- 10^floor(log10(max(conv.crit)))
    }
    else if (!is.null(next.tol.target) && max(conv.crit) <
             next.tol.target) {
      report.backfit.progress(verbose.lvl, next.tol.target)
      next.tol.target <- next.tol.target/10
    }

    ################################################ EDITED FOR ELBO-TIME TRACKING
    time.vec <- c(time.vec, Sys.time())
    elbo.vec <- c(elbo.vec, get.obj(flash))
    ################################################


  }
  if (method == "parallel") {
    parallel::stopCluster(cl)
  }
  if (iter == maxiter) {
    report.maxiter.reached(verbose.lvl)
  }
  if (get.obj(flash) > old.obj) {
    report.backfit.complete(verbose.lvl, get.obj(flash))
  }
  announce.wrapup(verbose.lvl)
  flash <- wrapup.flash(flash, output.lvl = 3L)
  report.completion(verbose.lvl)

  ################################################ EDITED FOR ELBO-TIME TRACKING
  time.vec <- as.numeric(time.vec[-1] - time.vec[1])

  out.list <- list(flash=flash, time.vec=time.vec, elbo.vec=elbo.vec)
  return(out.list)
  ################################################

  #return(flash)
}



# compile functions for a modest speed boost
################################################################################################
library(compiler)


alt.flash <- cmpfun(alt.flash)
ELBO <- cmpfun(ELBO)
flash.backfit.e <- cmpfun(flash.backfit.e)
G.PF <- cmpfun(G.PF)
INIT <- cmpfun(INIT)
MEAN.SOL.B <- cmpfun(MEAN.SOL.B)
SOL.B <- cmpfun(SOL.B)
SOL.G <- cmpfun(SOL.G)
SOL.Q <- cmpfun(SOL.Q)
SOL.REG.B <- cmpfun(SOL.REG.B)
SOL.TAU <- cmpfun(SOL.TAU)
