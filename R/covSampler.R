# this function uses gibbs sampling to estimate the posterior distribution
# of a sample's covariance matrix
# sources: https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution,
# Murphy, K. P. (2007). Conjugate bayesian analysis of the gaussian distribution (Tech. Rep.).
# University of British Columbia.

covSamp <- function(data, n.iter, n.burnin, thin, n.chains, pairwise, callback = function(){}, k0, df0){
  p <- ncol(data)

  c_post <- array(0, c(n.chains, n.iter, p, p))
  inds <- which(is.na(data), arr.ind = TRUE)
  dat_imp <- array(0, c(n.chains, n.iter, nrow(inds)))

  for (z in 1:n.chains) {

    if (pairwise) {
      dat_complete <- data
      # initial generation of complete data set with means as substitutes
      dat_complete[inds] <- colMeans(data, na.rm = TRUE)[inds[, 2]]
      # now the missing are being replaced in each iteration with draws from the conditional joints
      for (i in 1:n.iter) {
        pars <- preCompCovParams(dat_complete, k0, df0)
        cc <- sampleCov(pars)
        ms <- numeric(p)
        c_post[z, i, , ] <- cc
        # substitute missing values one by one, where each value is drawn conditional on the rest of the data
        cols <- unique(inds[, 2])
        for (ccc in cols) {
          rows <- inds[which(inds[, 2] == ccc), 1]
          mu1 <- ms[ccc]
          mu2 <- ms[-ccc]
          cc11 <- cc[ccc, ccc]
          cc21 <- cc[-ccc, ccc]
          cc12 <- cc[ccc, -ccc]
          cc22 <- cc[-ccc, -ccc]
          ccq <- cc11 - cc12 %*% solve(cc22) %*% cc21
          for (r in rows) {
            muq <- mu1 + cc12 %*% solve(cc22) %*% (as.numeric(dat_complete[r, -ccc]) - mu2)
            dat_complete[r, ccc] <- rnorm(1, muq, sqrt(ccq))
          }
        }
        dat_imp[z, i, ] <- dat_complete[inds]
        callback()
      }

    } else {
      pars <- preCompCovParams(data, k0, df0)
      for (i in 1:n.iter){
        c_post[z, i, , ] <- sampleCov(pars) # sample from inverse Wishart
        callback()
      }
    }
  }

  c_post_burned <- c_post[, (n.burnin + 1):n.iter, , , drop = FALSE]
  c_post_out <- c_post_burned[, seq(1, dim(c_post_burned)[2], thin), , , drop = FALSE]

  dat_imp_burned <- dat_imp[, (n.burnin + 1):n.iter, , drop = FALSE]
  dat_out <- dat_imp_burned[, seq(1, dim(dat_imp_burned)[2], thin), , drop = FALSE]


  return(list(cov_mat = c_post_out, dat_mis_samp_cov = dat_out))
}

preCompCovParams <- function(data, k0, df0) {
  n <- nrow(data)
  p <- ncol(data)
  # posterior covariance matrix ---------------------------------------------------
  k0 <- k0
  if (is.null(df0)) df0 <- p
  t <- diag(p)
  T0 <- diag(k0, nrow = p, ncol = p) # matrix inversion of diagonal matrix
  mu0 <- rep(0, p) # prior means
  vn <- df0 + n

  ym <- .colMeans(data, n, p)
  S <- cov(sweep(data, 2L, ym, `-`)) * (n - 1)

  Tn <- T0 + S + (k0 * n / (k0 + n)) * (ym - mu0) %*% t(ym - mu0)
  # drawing samples from posterior:
  Tn <- chol(chol2inv(chol(Tn)))
  dfChisq <- vn:(vn-p+1)
  utz <- upper.tri(matrix(0, p, p))

  return(list(vn = vn, Tn = Tn, p = p, dfChisq = dfChisq, utz = utz))
}

sampleCov <- function(par_list) {
  # sample from inverse Wishart
  cc <- rinvwishart2(par_list$vn, par_list$Tn, par_list$p, par_list$dfChisq, par_list$utz)
  return(cc)
}

# ------- customized covariance matrix sampling with cholesky decomposition -----------
rinvwishart2 <- function(nu, S, k = nrow(S), dfChisq = nu:(nu - k + 1),
                         utz = upper.tri(matrix(0, k, k))) {

  Z <- matrix(0, k, k)
  x <- rchisq(k, dfChisq)
  x[x == 0] <- 1e-100
  diag(Z) <- sqrt(x)
  if (k > 1) {
    Z[utz] <- rnorm(k * {k - 1} / 2)
  }
  return(chol2inv(Z %*% S))
}
