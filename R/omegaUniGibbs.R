# this function uses gibbs sampling to estimate the loadings and error variances
# of a cfa one factor model
# it returns the posterior distribution sample of omegas calculated from those parameters
# source: Lee, S.-Y. (2007). Structural equation modeling: A bayesian approach(Vol. 711). JohnWiley & Sons.
# p. 81 ff.
omegaSampler <- function(data, n.iter, n.burnin, thin, n.chains, pairwise, callback = function(){}){

  n <- nrow(data)
  p <- ncol(data)

  omm <- matrix(0, n.chains, n.iter)
  lll <- array(0, c(n.chains, n.iter, p))
  ppp <- array(0, c(n.chains, n.iter, p))

  inds <- which(is.na(data), arr.ind = TRUE)
  dat_imp <- array(0, c(n.chains, n.iter, nrow(inds)))

  # hyperparameters
  # prior multiplier for loadings variance, prior shape and rate for residuals, prior loadings,
  # prior scaling for cov matrix of factor scores, prior df for cov matrix of factor scores
  pars <- list(H0k = 1, a0k = 2, b0k = 1, l0k = numeric(p), R0 = p, p0 = p + 2)

  for (z in 1:n.chains) {
    # draw starting values for sampling from prior distributions:
    ss <- drawStart(n, p, pars)
    wi <- ss$wi
    phi <- ss$phi

    if (pairwise) { # missing data
      dat_complete <- data
      dat_complete[inds] <- colMeans(data, na.rm = TRUE)[inds[, 2]]
      ms <- rep(0, p)

      for (i in 1:n.iter) {
        out <- sampleFMParams(wi, dat_complete, phi, pars)
        wi <- out$wi
        phi <- out$phi
        cc <- out$cc
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
          ccq <- cc11 - cc12 %*% try(solve(cc22)) %*% cc21
          for (r in rows) {
            muq <- mu1 + cc12 %*% try(solve(cc22)) %*% (as.numeric(dat_complete[r, -ccc]) - mu2)
            dat_complete[r, ccc] <- rnorm(1, muq, sqrt(ccq))
          }
        }
        omm[z, i] <- omegaBasic(out$lambda, out$psi)
        dat_imp[z, i, ] <- dat_complete[inds]
        lll[z, i, ] <- out$lambda
        ppp[z, i, ] <- out$psi
        callback()
      }

    } else { # no missing data

      for (i in 1:n.iter){
        oo <- sampleFMParams(wi, data, phi, pars)
        omm[z, i] <- omegaBasic(oo$lambda, oo$psi)
        lll[z, i, ] <- oo$lambda
        ppp[z, i, ] <- oo$psi

        wi <- oo$wi
        phi <- oo$phi
        callback()
      }
    }
  }

  omm_burned <- omm[, (n.burnin + 1):n.iter, drop = FALSE]
  omm_out <- omm_burned[, seq(1, dim(omm_burned)[2], thin), drop = FALSE]

  lll_burned <- lll[, (n.burnin + 1):n.iter, , drop = FALSE]
  ppp_burned <- ppp[, (n.burnin + 1):n.iter, , drop = FALSE]
  lll_out <- lll_burned[, seq(1, dim(lll_burned)[2], thin), , drop = FALSE]
  ppp_out <- ppp_burned[, seq(1, dim(ppp_burned)[2], thin), , drop = FALSE]

  dat_imp_burned <- dat_imp[, (n.burnin + 1):n.iter, , drop = FALSE]
  dat_out <- dat_imp_burned[, seq(1, dim(dat_imp_burned)[2], thin), , drop = FALSE]


  return(list(omega = coda::mcmc(omm_out), lambda = coda::mcmc(lll_out), psi = coda::mcmc(ppp_out),
              dat_mis_samp_fm = coda::mcmc(dat_out)))
}



sampleFMParams <- function(wi, data, phi, pars) {
  n <- nrow(data)
  p <- ncol(data)

  H0k <- pars$H0k # prior multiplier for lambdas variance
  R0 <- pars$R0 # prior shape for wishart distribution for variance of factor scores (wi)
  p0 <- pars$p0 # prior df for wishart distribution for variance of factor scores (wi)

  l0k <- pars$l0k # prior lambdas
  a0k <- pars$a0k # prior shape parameter for gamma function for psis
  b0k <- pars$b0k # prior rate parameter for gamma for psi

  # hyperparameters for posteriors
  Ak <- (1 / H0k + c(t(wi) %*% wi))^-1
  ak <- Ak * ((1 / H0k) * l0k + t(wi) %*% data)
  bekk <- b0k + 0.5 * (t(data) %*% data - (t(ak) * (1 / Ak)) %*% ak
                       + (l0k * (1 / H0k)) %*% t(l0k))
  bek <- diag(bekk)

  #  sample psi and lambda
  invpsi <- rgamma(p, n / 2 + a0k, bek)
  invPsi <- diag(invpsi)
  psi <- 1 / invpsi
  lambda <- rnorm(p, ak * sqrt(as.vector(phi)), sqrt(psi * Ak))

  if (mean(lambda) < 0) {# solve label switching problem
    lambda <- -lambda
  }
  invphi <- 1 / phi
  # sample wi posterior:
  m <- solve(invphi + t(lambda) %*% invPsi %*% lambda) %*% t(lambda) %*% invPsi %*% t(data)
  V <- solve(invphi + t(lambda) %*% invPsi %*% lambda)
  wi <- rnorm(n, m, sqrt(V))
  # set factor variance to 1 to identify the model
  wi <- wi / sd(wi)

  # sample phi:
  phi <- LaplacesDemon::rinvwishart(nu = n + p0, S = t(wi) %*% (wi) + R0)
  invphi <- 1 / phi

  cc <- lambda %*% phi %*% t(lambda) + diag(psi) # phi = 1 is bad!

  return(list(psi = psi, lambda = lambda, phi = phi, wi = wi, cc = cc))
}


drawStart <- function(n, p, pars) {

  invpsi <- rgamma(p, pars$a0k, pars$b0k)
  psi <- 1 / invpsi

  lambda <- rnorm(p, pars$l0k, sqrt(psi * pars$H0k))

  phi <- LaplacesDemon::rinvwishart(nu = pars$p0, S = pars$R0)

  wi <- rnorm(n, 0, sqrt(phi))
  wi <- wi/sd(wi) # fix variance to 1

  return(list(psi = psi, lambda = lambda, wi = wi, phi = phi))
}
