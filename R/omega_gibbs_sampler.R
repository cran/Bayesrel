# this function uses gibbs sampling to estimate the loadings and error variances
# of a cfa one factor model
# it returns the posterior distribution sample of omegas calculated from those parameters
# source: Lee, S.-Y. (2007). Structural equation modeling: A bayesian approach(Vol. 711). JohnWiley & Sons.
# p. 81 ff.
omegaSampler <- function(data, n.iter, n.burnin, pairwise, thin = 1, n.chains = 1){

  n <- nrow(data)
  p <- ncol(data)
  H0k <- 1 # prior multiplier for lambdas variance

  l0k <- rep(0, p) # prior lambdas
  a0k <- 2 # prior shape parameter for gamma function for psis
  b0k <- 1 # prior rate parameter for gamma for psi


  R0 <- p # prior shape for wishart distribution for variance of factor scores (wi)
  p0 <- p+2 # prior df for wishart distribution for variance of factor scores (wi)
  # this lets the factor variance be approx 1

  omega <- matrix(0, (n.iter-n.burnin), n.chains)

  for (z in 1:n.chains) {
    # draw starting values for sampling from prior distributions:
    invpsi <- rgamma(p, a0k, b0k)
    psi <- 1/invpsi
    invPsi <- diag(invpsi)

    lambda <- rnorm(p, l0k, sqrt(psi*H0k))

    # phi <- 1
    phi <- LaplacesDemon::rinvwishart(nu = p0, S = R0)
    invphi <- 1/phi

    wi <- rnorm(n, 0, sqrt(phi))
    wi <- wi/sd(wi) # fix variance to 1

    # prepare matrices for saving lambda and psi and omega:
    La <- matrix(0, n.iter, p)
    Psi <- matrix(0, n.iter, p)
    oms <- numeric(n.iter)
    ph <- numeric(n.iter)


    if (pairwise) { # missing data
      inds <- which(is.na(data), arr.ind = T)
      dat_complete <- data
      dat_complete[inds] <- colMeans(data, na.rm = T)[inds[, 2]]

      for (i in 1:n.iter){

        # hyperparameters for posteriors
        Ak <- (1/H0k + c(t(wi) %*% wi))^-1
        ak <- Ak * ((1/H0k) * l0k + t(wi) %*% dat_complete)
        bekk <- b0k + 0.5 * (t(dat_complete) %*% dat_complete - (t(ak) * (1/Ak)) %*% ak
                             + (l0k * (1/H0k)) %*% t(l0k))
        bek <- diag(bekk)

        #  sample psi and lambda
        invpsi <- rgamma(p, n/2 + a0k, bek)
        invPsi <- diag(invpsi)
        psi <- 1/invpsi
        lambda <- rnorm(p, ak * sqrt(as.vector(phi)), sqrt(psi * Ak))

        if (mean(lambda) < 0) {# solve label switching problem
          lambda <- -lambda
        }
        # sample wi posterior:
        m <- solve(invphi + t(lambda) %*% invPsi %*% lambda) %*% t(lambda) %*% invPsi %*% t(dat_complete)
        V <- solve(invphi + t(lambda) %*% invPsi %*% lambda)
        wi <- rnorm(n, m, sqrt(V))
        # set factor variance to 1 to identify the model
        wi <- wi/sd(wi)

        # sample phi:
        phi <- LaplacesDemon::rinvwishart(nu = n + p0, S = t(wi) %*% (wi) + R0)
        invphi <- 1/phi

        # substitute missing values one by one, assuming a normal with mean zero and model implied variance
        # by conditioning on the remaining complete data:

        cc <- lambda %*% phi %*% t(lambda) + diag(psi)

        # ms <- MASS::mvrnorm(1, numeric(p), cc)
        ms <- numeric(p) # not sure if this is the right way, it works however
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

        oms[i] <- omegaBasic(lambda, psi)

        La[i, ] <- lambda
        Psi[i, ] <- psi

      }

    } else { # no missing data

      for (i in 1:n.iter){

        # hyperparameters for posteriors
        Ak <- (1/H0k + c(t(wi) %*% wi))^-1
        ak <- Ak * ((1/H0k) * l0k + t(wi) %*% data)
        bekk <- b0k + 0.5 * (t(data) %*% data - (t(ak) * (1/Ak)) %*% ak
                             + (l0k * (1/H0k)) %*% t(l0k))
        bek <- diag(bekk)

        #  sample psi and lambda
        invpsi <- rgamma(p, n/2 + a0k, bek)
        invPsi <- diag(invpsi)
        psi <- 1/invpsi
        lambda <- rnorm(p, ak * sqrt(as.vector(phi)), sqrt(psi * Ak))

        if (mean(lambda) < 0) # solve label switching problem
          lambda <- -lambda

        # sample wi posterior:
        m <- solve(invphi + t(lambda) %*% invPsi %*% lambda) %*% t(lambda) %*% invPsi %*% t(data)
        V <- solve(invphi + t(lambda) %*% invPsi %*% lambda)
        wi <- rnorm(n, m, sqrt(V))
        # set factor variance to 1 to identify the model
        wi <- wi/sd(wi)

        # sample phi:
        phi <- LaplacesDemon::rinvwishart(nu = n + p0, S = t(wi) %*% (wi) + R0)
        invphi <- 1/phi

        oms[i] <- omegaBasic(lambda, psi)

        La[i, ] <- lambda
        Psi[i, ] <- psi

      }
    }

    # n.burnin
    oms <- oms[(n.burnin + 1):n.iter]
    Psi <- Psi[(n.burnin + 1):n.iter, ]
    La <- La[(n.burnin + 1):n.iter, ]

    omega[, z] <- oms
  }
  omega <- as.vector(omega)
  omega <- omega[seq(1, length(omega), thin)]

  return(list(omega = coda::mcmc(omega), lambda = coda::mcmc(La), psi = coda::mcmc(Psi)
  ))
}
