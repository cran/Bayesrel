# this function samples priors for the estimates and the number of indicators

priorSampUni <- function(p, estimate, n.samp = 2e3, k0, df0, a0, b0, m0){
  group <- c("alpha", "lambda2", "lambda4", "lambda6", "glb")
  if (estimate %in% group) {
    v0 <- df0
    k0 <- k0
    t <- diag(p)
    T0 <- solve(t / k0)
    m <- array(0, c(n.samp, p, p))
    for (i in seq_len(n.samp)){
      m[i, , ] <- LaplacesDemon::rinvwishart(v0, T0)
    }
  }
  out <- list()
  if (estimate == "alpha") {
    prioralpha <- apply(m, MARGIN = 1, applyalpha)
    out <- density(prioralpha, from = 0, to = 1, n = 512)
  }
  if (estimate == "lambda2") {
    priorlambda2 <- apply(m, MARGIN = 1, applylambda2)
    out <- density(priorlambda2, from = 0, to = 1, n = 512)
  }
  if (estimate == "lambda4") {
    priorlambda4 <- apply(m, MARGIN = 1, applylambda4NoCpp)
    out <- density(priorlambda4, from = 0, to = 1, n = 512)
  }
  if (estimate == "lambda6") {
    priorlambda6 <- apply(m, MARGIN = 1, applylambda6)
    out <- density(priorlambda6, from = 0, to = 1, n = 512)
  }
  if (estimate == "glb") {
    priorglb <- glbOnArrayCustom(m)
    out <- density(priorglb, from = 0, to = 1, n = 512)
  }
  if (estimate == "omega") {
    H0 <- 1 # prior multiplier matrix for lambdas variance
    l0k <- rep(m0, p) # prior lambdas
    a0k <- a0 # prior gamma function for psis
    b0k <- b0 # prior gamma for psi
    prioromega <- numeric(n.samp)
    for (i in 1:n.samp) {
      invpsi <- rgamma(p, a0k, b0k)
      psi <- 1 / invpsi
      lambda <- rnorm(p, l0k, sqrt(psi * H0))
      prioromega[i] <- omegaBasic(lambda, psi)
    }
    out <- density(prioromega, from = 0, to = 1, n = 512)
  }

  return(out)

}


omegasSecoPrior <- function(x, nsamp = 2e3) {

  k <- x$k
  ns <- x$n.factors
  a0 <- x$prior_params[["a0"]]
  b0 <- x$prior_params[["b0"]]
  l0 <- x$prior_params[["l0"]]
  A0 <- x$prior_params[["A0"]]
  c0 <- x$prior_params[["c0"]]
  d0 <- x$prior_params[["d0"]]
  beta0 <- x$prior_params[["beta0"]]
  B0 <- x$prior_params[["B0"]]
  p0 <- x$prior_params[["p0"]]
  R0 <- x$prior_params[["R0"]]
  modelfile <- x$model$file

  idex <- modelfile$idex
  imat <- modelfile$imat
  # ---- sampling start --------
  l0mat <- matrix(0, k, ns)
  l0mat[imat] <- l0
  beta0vec <- numeric(ns)
  beta0vec[1:ns] <- beta0


  omh_prior <- numeric(nsamp)
  omt_prior <- numeric(nsamp)

  for (i in 1:nsamp) {

    phiw <- 1 / rgamma(1, shape = p0 / 2, scale = 2 / R0)

    invpsi <- rgamma(k, a0, b0)
    psi <- 1 / invpsi
    lambda <- rnorm(sum(imat), l0mat[imat], sqrt(psi * A0))
    # structural parameters
    invpsiw <- rgamma(ns, c0, d0)
    psiw <- 1 / invpsiw
    beta <- rnorm(ns, beta0vec * sqrt(phiw), sqrt(psiw * B0))

    lmat <- l0mat
    lmat[imat] <- lambda
    lmat <- cbind(0, lmat)
    bmat <- matrix(0, ns + 1, ns + 1)
    bmat[2:(ns + 1), 1] <- beta
    om_prior <- omegasSeco(lmat, bmat, diag(psi), diag(c(1, psiw)))
    omh_prior[i] <- om_prior[1]
    omt_prior[i] <- om_prior[2]

  }

  return(list(omh_prior = omh_prior, omt_prior = omt_prior))

}


omegasBifPrior <- function(x, nsamp = 2e3) {

  k <- x$k
  ns <- x$n.factors
  a0 <- x$prior_params[["a0"]]
  b0 <- x$prior_params[["b0"]]
  l0 <- x$prior_params[["l0"]]
  A0 <- x$prior_params[["A0"]]
  beta0 <- x$prior_params[["beta0"]]
  B0 <- x$prior_params[["B0"]]
  p0 <- x$prior_params[["p0"]]
  R0 <- x$prior_params[["R0"]]
  modelfile <- x$model$file

  idex <- modelfile$idex
  imat <- modelfile$imat
  # ---- sampling start --------
  l0mat <- matrix(0, k, ns)
  l0mat[imat] <- l0
  beta0vec <- numeric(k)
  beta0vec[1:k] <- beta0


  omh_prior <- numeric(nsamp)
  omt_prior <- numeric(nsamp)

  for (i in 1:nsamp) {

    phi <- 1 / rgamma(1, shape = p0 / 2, scale = 2 / R0)

    invpsi <- rgamma(k, a0, b0)
    psi <- 1 / invpsi

    imatb <- cbind(TRUE, imat)
    lmat <- cbind(beta0vec, l0mat)
    m0 <- lmat[imatb]

    lambda <- rnorm(length(m0), m0 * sqrt(phi), sqrt(c(psi * B0, psi * A0)))
    lmat[imatb] <- lambda

    om_prior <- omegasBif(lmat[, -1], lmat[, 1], diag(psi))
    omh_prior[i] <- om_prior[1]
    omt_prior[i] <- om_prior[2]

  }

  return(list(omh_prior = omh_prior, omt_prior = omt_prior))
}


omegasCorrPrior <- function(x, nsamp = 2e3) {

  k <- x$k
  ns <- x$n.factors
  a0 <- x$prior_params[["a0"]]
  b0 <- x$prior_params[["b0"]]
  l0 <- x$prior_params[["l0"]]
  A0 <- x$prior_params[["A0"]]
  p0 <- x$prior_params[["p0"]]
  R0 <- x$prior_params[["R0"]]
  modelfile <- x$model$file

  idex <- modelfile$idex
  imat <- modelfile$imat
  # ---- sampling start --------
  l0mat <- matrix(0, k, ns)
  l0mat[imat] <- l0

  H0k <- rep(A0, ns)


  omt_prior <- numeric(nsamp)

  for (i in 1:nsamp) {

    phi <- LaplacesDemon::rinvwishart(nu = p0, S = R0)

    invpsi <- rgamma(k, a0, b0)
    psi <- 1 / invpsi
    lambda <- rnorm(sum(imat), l0mat[imat], sqrt(psi * A0))

    lmat <- l0mat
    lmat[imat] <- lambda

    om_prior <- omegaCorr(lmat, phi, diag(psi))
    omt_prior[i] <- om_prior

  }

  return(list(omt_prior = omt_prior))

}

