library(progress)

genX <- function(groups,
                 group_ratio,
                 signal_lvl,
                 n,
                 p,
                 F_dist,
                 F_sd,
                 noise_sd,
                 seed=1
){

  K <- ncol(groups)

  L_true <- groups[ rep(1:length(group_ratio), times = n * group_ratio), ]
  for (k in 1:K){ L_true[, k] <- L_true[, k] * signal_lvl[k] }

  set.seed(seed)
  if (F_dist=='normal') {
    F_vals <- rnorm(p * K, mean = 0, sd = F_sd)
  } else if(F_dist=='laplace'){
    F_vals <- distr::r(distr::DExp(rate = sqrt(2)/F_sd))(p * K)
  }

  F_true <- matrix(F_vals, nrow = p, ncol = K)
  X <- tcrossprod(L_true, F_true) +
    matrix(rnorm(n * p, mean = 0, sd = noise_sd), nrow = n, ncol = p)

  out.list <- list(X = X, L_true = L_true, F_true = F_true)
  return(out.list)
}



# input: X.list, (ebcd or flash) fit$L.pm
# output: cs (cosine similarity),
#         scaled and column-sorted L estimate (that achieves maximum cs)

eval.fit <- function(X.list, L.pm){

  if (is.null(L.pm)){
    out.list <- list(cs = NA, L_est = NA)
  } else if (ncol(L.pm) != ncol(X.list$L_true)){
    out.list <- list(cs = NA, L_est = NA)
  } else{
  L_est <- scale(L.pm, center = FALSE, scale = TRUE)

  matchings <- combinat::permn(ncol(X.list$L_true))
  L_true_colnorm <- sqrt(colSums(X.list$L_true^2))
  L_est_colnorm <- sqrt(colSums(L_est^2))

  cs <- rep(0, length(matchings))
  for (m in 1:length(matchings)) {
    cs[m] <- mean(abs(colSums(L_est[, matchings[[m]]] * X.list$L_true)) / (
      L_est_colnorm[matchings[[m]]] * L_true_colnorm))
  }
  L_est <- L_est[, matchings[[which.max(cs)]]]
  cs <- max(cs)

  out.list <- list(cs = cs, L_est = L_est)
  }
  return(out.list)
}


# iterations
n.cs <- function(genX.iter, method, ebnm.fn, nseed = 50){

  t0 <- Sys.time()

  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)",
                         total = nseed)
  cs.all <- rep(0, nseed)
  for (seed in 1:nseed){

    pb$tick()
    X.list <- genX.iter(seed = seed)
    if(method=='ebcd'){
      fit <- ebcd(data = X.list$X,
                  ebnm.fn = ebnm.fn,
                  greedy.Kmax = ncol(X.list$L_true), # use true K
                  verbose = 0)
    }else if(method=='flash'){
      fit <- flash(data = X.list$X,
                   ebnm.fn = ebnm.fn,
                   greedy.Kmax = ncol(X.list$L_true), # use true K
                   backfit = TRUE,
                   verbose = 0)
    }else if(method=='covflash'){
      fit <- flash(data = tcrossprod(X.list$X),
                   ebnm.fn = ebnm.fn,
                   greedy.Kmax = ncol(X.list$L_true), # use true K
                   backfit = TRUE,
                   verbose = 0)
    }

    cs.all[seed] <- eval.fit(X.list, fit$L.pm)$cs
  }

  t1 <- Sys.time()

  out.list <- list(cs = cs.all,
                   time = t1 - t0)
  return(out.list)
}



