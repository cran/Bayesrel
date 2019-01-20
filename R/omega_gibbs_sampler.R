# this function uses gibbs sampling to estimate the loadings and error variances
# of a cfa one factor model
# it returns the posterior distribution sample of omegas calculated from those parameters
# source: Lee, S.-Y. (2007). Structural equation modeling: A bayesian approach(Vol. 711). JohnWiley & Sons.
# p. 81 ff.
omegaSampler <- function(data, n.iter = 2e3, n.burnin = 50){
  n <- nrow(data)
  p <- ncol(data)
  H0 <- 2.5 # prior multiplier matrix for lambdas variance

  l0k <- rep(0, p) # prior lambdas
  a0k <- 1 # prior gamma function for psis
  b0k <- .05 # prior gamma for psi

  # draw starting values for sampling from prior distributions:
  invpsi <- rgamma(p, a0k, b0k)
  invPsi <- diag(invpsi)
  psi <- 1/invpsi
  lambda <- rnorm(p, l0k, sqrt(psi * H0))

  Phi <- 1
  invPhi <- 1/Phi
  wi <- rnorm(n, 0, sqrt(Phi))
  wi <- wi/sd(wi) # fix variance to 1
  # prepare matrices for saving lambda and psi:
  La <- matrix(0, n.iter, p)
  Ps <- matrix(0, n.iter, p)

  for (i in 1:n.iter){
    # hyperparameters for posteriors
    Ak <- as.vector(1/((1/H0) + t(wi) %*% wi))
    ak <- as.vector(Ak * ((1/H0) * l0k + t(wi) %*% data))
    bekk <- b0k + 0.5 * (t(data) %*% data - ak %*% t(ak)* (1/Ak)
                         + l0k %*% t(l0k) * (1/H0))
    bek <- diag(bekk)
    #  sample psi and lambda
    invpsi <- rgamma(p, n/2 + a0k, bek)
    invPsi <- diag(invpsi)
    psi <- 1/invpsi
    lambda <- rnorm(p, ak, sqrt(psi * Ak))

    # sample wi posterior:
    m <- solve(invPhi + t(lambda) %*% invPsi %*% lambda) %*% t(lambda) %*% invPsi %*% t(data)
    V <- solve(invPhi + t(lambda) %*% invPsi %*% lambda)
    wi <- rnorm(n, m, sqrt(V))
    wi <- wi/sd(wi)
    # we dont sample phi

    La[i, ] <- lambda
    Ps[i, ] <- psi
  }
  # n.burnin
  La <- La[(n.burnin + 1):n.iter, ]
  Ps <- Ps[(n.burnin + 1):n.iter, ]

  gibbs.o.obj <- numeric(nrow(La))
  for (i in 1:(nrow(La))){
    gibbs.o.obj[i] <- sum(La[i,])^2 / (sum(La[i,])^2 + sum(Ps[i,]))
  }
  return(list(omega = gibbs.o.obj, lambda = La, psi = Ps))
}

