# this function samples priors for the estimates and the number of indicators

priorSampUni <- function(p, estimates, n.samp = 2e3){
  if ("alpha" %in% estimates || "lambda2" %in% estimates || "lambda4" %in% estimates || "lambda6" %in% estimates ||
      "glb" %in% estimates){
    v0 <- p
    k0 <- 1e-10
    t <- diag(p)
    T0 <- solve(t / k0)
    m <- array(0, c(n.samp, p, p))
    for (i in 1:n.samp){
      m[i, , ] <- LaplacesDemon::rinvwishart(v0, T0)
    }
  }
  out <- list()
  if ("alpha" %in% estimates) {
    priora <- apply(m, MARGIN = 1, applyalpha)
    out$prioralpha <- quantiles(priora[priora >= 0])
  }
  if ("lambda2" %in% estimates) {
    priorlambda2 <- apply(m, MARGIN = 1, applylambda2)
    out$priorlambda2 <- quantiles(priorlambda2[priorlambda2 >= 0])
  }
  if ("lambda4" %in% estimates) {
    priorlambda4 <- apply(m, MARGIN = 1, applylambda4NoCpp)
    out$priorlambda4 <- quantiles(priorlambda4[priorlambda4 >= 0])
  }
  if ("lambda6" %in% estimates) {
    priorlambda6 <- apply(m, MARGIN = 1, applylambda6)
    out$priorlambda6 <- quantiles(priorlambda6[priorlambda6 >= 0])
  }
  if ("glb" %in% estimates) {
    priorglb <- glbOnArrayCustom(m)
    out$priorglb <- quantiles(priorglb[priorglb >= 0])
  }
  if ("omega" %in% estimates) {
    H0 <- 1 # prior multiplier matrix for lambdas variance
    l0k <- rep(0, p) # prior lambdas
    a0k <- 1 # prior gamma function for psis
    b0k <- 2 # prior gamma for psi
    prioromega <- numeric(n.samp)
    for (i in 1:n.samp) {
      invpsi <- rgamma(p, a0k, b0k)
      psi <- 1 / invpsi
      lambda <- rnorm(p, l0k, sqrt(psi * H0))
      prioromega[i] <- omegaBasic(lambda, psi)
    }
    out$prioromega <- quantiles(prioromega[prioromega >= 0])
  }

  return(out)

}


omegasPrior <- function(k, ns, nsamp = 2e3) {

  # index matrix for lambdas
  idex <- matrix(seq(1:k), ns, k / ns, byrow = TRUE)
  imat <- matrix(FALSE, k, ns)
  for (i in 1:ns) {
    imat[idex[i, ], i] <- TRUE
  }

  H0k <- rep(1, ns) # prior multiplier matrix for lambdas variance (and covariance)
  l0k <- matrix(0, k, ns) # prior lambdas
  a0k <- 2 # prior shape parameter for gamma function for psis
  b0k <- 1 # prior rate parameter for gamma for psi

  # -------------- structural equation -----------
  H0kw <- 2.5
  beta0k <- numeric(ns)
  a0kw <- 2
  b0kw <- 1

  R0w <- diag(rep(1 / (k), ns + 1))
  p0w <- ns^2

  # ---- sampling start --------

  pars <- list(H0k = H0k, a0k = a0k, b0k = b0k, l0k = l0k,
               H0kw = H0kw, a0kw = a0kw, b0kw = b0kw, beta0k = beta0k,
               R0w = R0w, p0w = p0w)
  omh_prior <- numeric(nsamp)
  omt_prior <- numeric(nsamp)

  for (i in 1:nsamp) {
    invpsi <- rgamma(k, pars$a0k, pars$b0k)
    psi <- 1 / invpsi
    lambda <- rnorm(k, pars$l0k[imat], sqrt(psi * rep(pars$H0k, each = k / ns)))
    # structural parameters
    invpsiw <- rgamma(ns, pars$a0kw, pars$b0kw)
    psiw <- 1 / invpsiw
    beta <- rnorm(ns, pars$beta0k, sqrt(psiw*pars$H0kw))

    lmat <- pars$l0k
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
