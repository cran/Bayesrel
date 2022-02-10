# this function samples priors for the estimates and the number of indicators

priorSampUni <- function(p, estimates, n.samp = 2e3, k0, df0, a0, b0, m0){
  if ("alpha" %in% estimates || "lambda2" %in% estimates || "lambda4" %in% estimates || "lambda6" %in% estimates ||
      "glb" %in% estimates){
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
  if ("alpha" %in% estimates) {
    prioralpha <- apply(m, MARGIN = 1, applyalpha)
    out <- density(prioralpha, from = 0, to = 1, n = 512)
  }
  if ("lambda2" %in% estimates) {
    priorlambda2 <- apply(m, MARGIN = 1, applylambda2)
    out <- density(priorlambda2, from = 0, to = 1, n = 512)
  }
  if ("lambda4" %in% estimates) {
    priorlambda4 <- apply(m, MARGIN = 1, applylambda4NoCpp)
    out <- density(priorlambda4, from = 0, to = 1, n = 512)
  }
  if ("lambda6" %in% estimates) {
    priorlambda6 <- apply(m, MARGIN = 1, applylambda6)
    out <- density(priorlambda6, from = 0, to = 1, n = 512)
  }
  if ("glb" %in% estimates) {
    priorglb <- glbOnArrayCustom(m)
    out <- density(priorglb, from = 0, to = 1, n = 512)
  }
  if ("omega" %in% estimates) {
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


omegasPrior <- function(k, ns, nsamp = 2e3,
                        a0, b0, l0, A0, c0, d0, beta0, B0, p0, R0) {

  # index matrix for lambdas
  idex <- matrix(seq(1:k), ns, k / ns, byrow = TRUE)
  imat <- matrix(FALSE, k, ns)
  for (i in 1:ns) {
    imat[idex[i, ], i] <- TRUE
  }

  # ---- sampling start --------
  l0mat <- matrix(0, k, ns)
  l0mat[imat] <- l0
  beta0vec <- numeric(ns)
  beta0vec[1:ns] <- beta0

  pars <- list(H0k = rep(A0, ns), a0k = a0, b0k = b0, l0k = l0mat,
               H0kw = B0, a0kw = c0, b0kw = d0, beta0k = beta0vec,
               p0w = p0, R0winv = R0)
  omh_prior <- numeric(nsamp)
  omt_prior <- numeric(nsamp)

  for (i in 1:nsamp) {

    phiw <- 1 / rgamma(1, shape = pars$p0w / 2, scale = 2 / pars$R0winv)

    invpsi <- rgamma(k, pars$a0k, pars$b0k)
    psi <- 1 / invpsi
    lambda <- rnorm(k, pars$l0k[imat], sqrt(psi * rep(pars$H0k, each = k / ns)))
    # structural parameters
    invpsiw <- rgamma(ns, pars$a0kw, pars$b0kw)
    psiw <- 1 / invpsiw
    beta <- rnorm(ns, pars$beta0k * sqrt(phiw), sqrt(psiw * pars$H0kw))

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
