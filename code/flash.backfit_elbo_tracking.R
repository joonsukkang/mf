function (flash, kset = NULL, method = c("extrapolate", "sequential",
                                         "dropout", "random", "parallel"), warmstart = TRUE, conv.crit.fn = calc.obj.diff,
          tol = set.default.tol(flash), maxiter = 500, verbose.lvl = get.verbose.lvl(flash))
{

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
