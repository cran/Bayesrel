# forces a quadratic matrix to be symmetrical
make_symmetric <- function(a, lower.tri = TRUE) {
  if (lower.tri){
    ind <- upper.tri(a)
    a[ind] <- t(a)[ind]
  } else {
    ind <- lower.tri(a)
    a[ind] <- t(a)[ind]
  }
  return(a)
}

# computes alpha analytical interval with given bounds
# source:
# Bonett, D. G., & Wright, T. A. (2015). Cronbach’s alpha reliability:
# Interval estimation, hypothesis testing, and sample size planning. Journal of Organizational Behavior,36(1), 3–15.
#
ciAlpha <- function(palpha, n, V) {
  p <- ncol(V)
  z <- qnorm(1 - palpha/2)
  b <- log(n/ (n - 1))
  j <- matrix(rep(1, p), nrow = p, ncol = 1)
  a0 <- t(j) %*% V %*% j
  t1 <- sum(diag(V))
  t2 <- sum(diag(V %*% V))
  a1 <- a0 ^ 3
  a2 <- a0 * (t2 + t1^2)
  a3 <- 2 * t1 * t(j) %*% (V %*% V) %*% j
  a4 <- 2 * p^2 / (a1 * (p - 1)^2)
  r <- (p/ (p-1)) * (1 - t1 / a0)
  var <- a4 * (a2 - a3) / (n - 3)
  ll <- 1 - exp(log(1 - r) - b + z * sqrt(var / (1 - r)^2))
  ul <- 1 - exp(log(1 - r) - b - z * sqrt(var / (1 - r)^2))
  out <- c(ll, ul)
  return(out)
}

# does quantile averaging and returns 2000 datapoints
quantiles <- function(samp, length_out = 2e3){
  q <- quantile(samp, probs = seq(0, 1, length.out = length_out))
  return(q)
}


se <- function(x) {
  b <- length(x)
  se <- sqrt(1 / (b-1) * sum((x - mean(x))^2))
  se
}




# calculate the kolomogorov smirnov distances between some samples and the original sample
KSTestStatistic <- function(x, y) {
  t <- stats::ks.test(x, y)
  return(t$statistic)
}

# calculate the kublack leibler distance between two samples
KLDStatistic <- function(x, y) {
  # transform the samples to PDFs:
  xdf <- getApproxDensity(x)
  ydf <- getApproxDensity(y)

  xx <- seq(0, 1, length.out = 1e3)
  t <- LaplacesDemon::KLD(xdf(xx), ydf(xx))
  return(t$sum.KLD.py.px)
}

hpdHelp <- function(x) {
  x <- coda::as.mcmc(x)
  return(coda::HPDinterval(x))
}


# create covariance matrix
createUnidimCovMat <- function(avg, p) {
  mean_cor <- 1
  counter <- 1
  while (mean_cor < (avg - .001) || mean_cor > (avg + .001)) {
    mlam <- avg * 3 + .02
    vlam <- avg * 2
    lam_true <- abs(rnorm(p, mlam, vlam))
    psi_true <- 1/rgamma(p, 10, 10)
    loading <- matrix(lam_true, nrow = p)
    psi_m <- diag(1, nrow = p)
    diag(psi_m) <- psi_true
    tmp_cov <- make_symmetric(loading %*% 1 %*% t(loading) + psi_m)
    cormat <- cov2cor(tmp_cov)
    mean_cor <- (sum(cormat) - p) / (p*p - p)
    counter <- counter + 1
    if (counter == 1e4)
      return(print("solution has not been found"))
  }
  return(tmp_cov)
}

trySmc <- function(M) {
  return(try(1 - 1 / diag(solve(cov2cor(M))), silent = TRUE))
}

tryPsd <- function(M) {
  R <- cov2cor(M)
  return(try(eigen(R, only.values = TRUE)$values, silent = TRUE))
}


getApproxDensity <- function(x) {
  d <- density(x, n = 2^12)
  f <- approxfun(d$x, d$y, yleft = 0, yright = 0)
  c <- integrate(f, 0, 1)$value
  return(
    function(x) {
      return(f(x) / c)
    }
  )
}

checkMcmcSettings <- function(n.iter, n.burnin, n.chains, thin) {

  if (n.chains < 2) stop("The number of chains needs to be at least 2.")
  if (n.iter <= n.burnin) stop("Too few samples remain after burnin.")
  if (length(seq(1, n.iter - n.burnin, thin)) < 100) stop("Too few posterior samples would result with the current MCMC
                                                          settings. Increase the iterations.")
}
