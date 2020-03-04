# this function uses gibbs sampling to estimate the posterior distribution
# of a sample's covariance matrix
# sources: https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution,
# Murphy, K. P. (2007). Conjugate bayesian analysis of the gaussian distribution(Tech. Rep.). University of British Columbia.

covSamp <- function(data, n.iter, n.burnin, pairwise){
  n <- nrow(data)
  p <- ncol(data)
  # posterior covariance matrix ---------------------------------------------------
  k0 <- 1e-10
  v0 <- p
  t <- diag(p)
  T0 <- diag(k0, nrow = p, ncol = p) # matrix inversion of diagonal matrix
  mu0 <- rep(0, p) # prior means
  kn <- k0 + n
  vn <- v0 + n
  c_post <- array(0, c(n.iter, p, p))

  if (pairwise) {
    inds <- which(is.na(data), arr.ind = T)
    dat_complete <- data
    # initial generation of complete data set with means as substitutes
    dat_complete[inds] <- colMeans(data, na.rm = T)[inds[, 2]]
    # now the missing are being replaced in each iteration with draws from the conditional joints
    for (i in 1:n.iter) {
      ym <- .colMeans(dat_complete, n, p)
      mun <- (k0 * mu0 + n * ym) / (k0 + n)
      S <- 0
      for (z in 1:n){
        S <- S + tcrossprod(dat_complete[z, ] - ym)
      }
      Tn <- T0 + S + (k0 * n / (k0 + n)) * (ym - mu0) %*% t(ym - mu0)
      # drawing samples from posterior:
      Tn <- chol(chol2inv(chol(Tn)))
      dfChisq <- vn:(vn-p+1)
      utz <- upper.tri(matrix(0, p, p))
      cc <- rinvwishart2(vn, Tn, p, dfChisq, utz) # sample from inverse Wishart
      # ms <- MASS::mvrnorm(1, mun, cc/kn)
      ms <- numeric(p)
      c_post[i, , ] <- cc
      # substitute missing values one by one, where each value is drawn conditional on the rest of the data
      rows <- unique(inds[, 1])
      for (r in rows) {
        cols <- inds[which(inds[, 1] == r), 2]
        for (b in cols) {
          mu1 <- ms[b]
          mu2 <- ms[-b]
          cc11 <- cc[b, b]
          cc21 <- cc[-b, b]
          cc12 <- cc[b, -b]
          cc22 <- cc[-b, -b]
          muq <- mu1 + cc12 %*% solve(cc22) %*% (as.numeric(dat_complete[r, -b]) - mu2)
          ccq <- cc11 - cc12 %*% solve(cc22) %*% cc21
          dat_complete[r, b] <- rnorm(1, muq, ccq)
        }
      }
    }

  } else {
    ym <- .colMeans(data, n, p)
    mun <- (k0 * mu0 + n * ym) / (k0 + n)
    S <- 0
    for (i in 1:n){
      S <- S + tcrossprod(data[i, ] - ym)
    }
    Tn <- T0 + S + (k0 * n / (k0 + n)) * (ym - mu0) %*% t(ym - mu0)
    # drawing samples from posterior:
    Tn <- chol(chol2inv(chol(Tn)))
    dfChisq <- vn:(vn-p+1)
    utz <- upper.tri(matrix(0, p, p))
    for (i in 1:n.iter){
      c_post[i, , ] <- rinvwishart2(vn, Tn, p, dfChisq, utz) # sample from inverse Wishart
    }
  }

  c_post <- c_post[(n.burnin + 1):n.iter, , ]

  return(c_post)
}

# ------- customized covariance matrix sampling with cholesky decomposition -----------
rinvwishart2 <- function(nu, S, k = nrow(S), dfChisq = nu:(nu-k+1), utz = upper.tri(matrix(0, k, k))) {

  # LaplacesDemon::rwishartc
  Z <- matrix(0, k, k)
  x <- rchisq(k, dfChisq)
  x[x == 0] <- 1e-100
  diag(Z) <- sqrt(x)
  if (k > 1) {
    # kseq <- 1:(k - 1)
    # Z[rep(k * kseq, kseq) + unlist(lapply(kseq, seq))] <- rnorm(k * {k - 1} / 2)
    # --end of copied code
    Z[utz] <- rnorm(k * {k - 1} / 2)
  }
  # LaplacesDemon::rinvwishart
  return(chol2inv(Z %*% S))
}


# covSamp_old <- function(data, n.iter, n.burnin){
#   n <- nrow(data)
#   p <- ncol(data)
#   # posterior covariance matrix ---------------------------------------------------
#   k0 <- 1e-10
#   v0 <- p
#   t <- diag(p)
#   T0 <- solve(t/k0) # inverse scale matrix, prior
#   mu0 <- rep(0, p) # prior means
#   ym <- apply(data, 2, mean)
#   # https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution, murphy 2007
#   mun <- (k0 * mu0 + n * ym) / (k0 + n)
#   kn <- k0 + n
#   vn <- v0 + n
#   S <- 0
#   for (i in 1:n){
#     M <- (data[i, ] - ym) %*% t(data[i, ] - ym)
#     S <- S + M
#   }
#   Tn <- T0 + S + (k0 * n / (k0 + n)) * (ym - mu0) %*% t(ym - mu0)
#   # drawing samples from posterior:
#   c_post <- array(0, c(n.iter, p, p))
#   for ( i in 1:n.iter){
#     c_post[i, , ] <- LaplacesDemon::rinvwishart(vn, Tn) # sample from inverse wishart
#   }
#   c_post <- c_post[(n.burnin + 1):n.iter, , ]
#
#   return(c_post)
# }
