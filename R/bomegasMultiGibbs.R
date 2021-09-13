


omegaMultiB <- function(data, ns, n.iter, n.burnin, n.chains, thin, model, pairwise, callback) {

  n <- nrow(data)
  k <- ncol(data)

  # ---- get the index matrix aka which items load on which factors ----
  mod_opts <- indexMatrix(model, k, ns, colnames(data))
  idex <- mod_opts$idex
  imat <- mod_opts$imat
  # ---- check missingsness ----
  inds <- which(is.na(data), arr.ind = TRUE)
  imputed <- array(0, c(n.chains, n.iter, nrow(inds)))

  # ---- sampling start --------

  pars <- list(H0k = rep(1, ns), a0k = 2, b0k = 1, l0k = matrix(0, k, ns),
               H0kw = 2.5, a0kw = 2, b0kw = 1, beta0k = numeric(ns),
               R0w = diag(rep(1 / k, ns + 1)), p0w = ns^2)

  omsh <- matrix(0, n.chains, n.iter)
  omst <- matrix(0, n.chains, n.iter)
  impl_covs <- array(0, c(n.chains, n.iter, k, k))

  for (ai in 1:n.chains) {
    # draw starting values
    starts <- drawStartMulti(n, k, ns, pars, imat)
    wi <- starts$wi
    phiw <- starts$phiw

    if (pairwise) { # missing data
      dat_filled <- data
      dat_filled[inds] <- colMeans(data, na.rm = TRUE)[inds[, 2]]
      ms <- rep(0, k)

      for (i in 1:n.iter) {
        params <- sampleSecoParams(dat_filled, pars, wi, phiw, ns, idex)
        wi <- params$wi
        phiw <- params$phiw
        # compute omega
        Lm <- cbind(0, params$lambda)
        Bm <- matrix(0, ns + 1, ns + 1)
        Bm[2:(ns + 1), 1] <- params$beta
        oms <- omegasSeco(Lm, Bm, diag(params$psi), diag(c(1, params$psiw)))

        omsh[ai, i] <- oms[1]
        omst[ai, i] <- oms[2]

        cc <- implCovMulti(Lm, Bm, theta = diag(params$psi), psi = diag(c(1, params$psiw)))
        impl_covs[ai, i, , ] <- cc

        # substitute missing values one by one, where each value is drawn conditional on the rest of the data
        # see https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
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
            muq <- mu1 + cc12 %*% try(solve(cc22)) %*% (as.numeric(dat_filled[r, -ccc]) - mu2)
            dat_filled[r, ccc] <- rnorm(1, muq, sqrt(ccq))
          }
        }
        imputed[ai, i, ] <- dat_filled[inds]

      }

    } else {
      for (i in 1:n.iter) {
        params <- sampleSecoParams(data, pars, wi, phiw, ns, idex)
        wi <- params$wi
        phiw <- params$phiw

        # compute omega
        Lm <- cbind(0, params$lambda)
        Bm <- matrix(0, ns + 1, ns + 1)
        Bm[2:(ns + 1), 1] <- params$beta
        oms <- omegasSeco(Lm, Bm, diag(params$psi), diag(c(1, params$psiw)))

        omsh[ai, i] <- oms[1]
        omst[ai, i] <- oms[2]

        impl_covs[ai, i, , ] <- implCovMulti(Lm, Bm, theta = diag(params$psi), psi = diag(c(1, params$psiw)))

      }
    }
  }

  # burnin
  omh_burn <- omsh[, (n.burnin + 1):n.iter, drop = F]
  omt_burn <- omst[, (n.burnin + 1):n.iter, drop = F]
  impl_covs_burn <- impl_covs[, (n.burnin + 1):n.iter, , , drop = F]

  # thinning
  omh_out <- omh_burn[, seq(1, dim(omh_burn)[2], thin), drop = F]
  omt_out <- omt_burn[, seq(1, dim(omt_burn)[2], thin), drop = F]
  impl_covs_out <- impl_covs_burn[, seq(1, dim(omt_burn)[2], thin), , , drop = F]

  return(list(omh = omh_out, omt = omt_out, impl_covs = impl_covs_out, imputed_values = imputed,
              modfile = mod_opts))

}


sampleSecoParams <- function(data, pars, wi, phiw, ns, idex) {

  n <- nrow(data)
  k <- ncol(data)
  # -------- measurement equation ------------
  H0k <- pars$H0k # prior multiplier matrix for lambdas variance (and covariance)
  l0k <- pars$l0k # prior lambdas
  a0k <- pars$a0k# prior shape parameter for gamma function for psis
  b0k <- pars$b0k  # prior rate parameter for gamma for psi

  # -------------- structural equation -----------
  H0kw <- pars$H0kw
  beta0k <- pars$beta0k
  a0kw <- pars$a0kw
  b0kw <- pars$b0kw

  R0w <- pars$R0w
  p0w <- pars$p0w

  ll <- matrix(0, k, ns)
  pp <- numeric(k)

  for (ii in 1:ns) {
    ids <- idex[[ii]]
    Ak <- solve(1 / H0k[ii] + t(wi[, ii + 1]) %*% wi[, ii + 1])
    ak <- Ak %*% (c(1 / H0k[ii]) %*% t(l0k[ids, ii]) + wi[, ii + 1] %*% data[, ids])
    bekk <- b0k + 0.5 * (t(data[, ids]) %*% data[, ids]
                         - t(ak) %*% solve(Ak) %*% ak
                         + (l0k[ids, ii] * 1 / H0k[ii]) %*% t(l0k[ids, ii]))
    bek <- diag(bekk)
    invpsi <- rgamma(length(ids), n / 2 + a0k, bek)
    psi <- 1 / invpsi
    lambda <- rnorm(length(ids), ak, sqrt(psi * as.vector(Ak)))

    if (mean(lambda) < 0) {# solve label switching problem
      lambda <- -lambda
    }
    ll[ids, ii] <- lambda
    pp[ids] <- psi
  }


  # ------- structural equation -----
  Akw <- 1 / (1 / H0kw + c(t(wi[, 1]) %*% wi[, 1]))
  akw <- Akw * (1 / H0kw * beta0k + t(wi[, 1]) %*% wi[, 2:(ns + 1)])
  bekkw <- b0kw + 0.5 * (t(wi[, 2:(ns + 1)]) %*% wi[, 2:(ns + 1)]
                         - ((t(akw) * (1 / Akw)) %*% akw)
                         + (beta0k * (1 / H0kw)) %*% t(beta0k))
  bekw <- diag(bekkw)

  invpsiw <- rgamma(ns, n / 2 + a0kw, bekw)
  psiw <- 1 / invpsiw
  beta <- rnorm(ns, akw * sqrt(diag(phiw)[1]), sqrt(psiw * Akw))

  if (mean(beta) < 0) {# solve label switching problem
    beta <- -beta
  }

  # in Lee it says to replace the usual inv Phi matrix when sampling the factor scores with a function
  # of the g-factor loadings and their residuals
  betaMat <- matrix(0, ns + 1, ns + 1)
  betaMat[2:(ns + 1), 1] <- beta
  ident <- diag(ns + 1)
  identInv <- solve(ident - betaMat)
  psid <- diag(c(1, psiw))
  sigW <- identInv %*% psid %*% t(identInv)
  invsigW <- solve(sigW)

  lll <- cbind(0, ll)
  impM <- t(lll) %*% solve(diag(pp)) %*% lll

  mw <- solve(invsigW + impM) %*% t(lll) %*% solve(diag(pp)) %*% t(data)
  Vw <- solve(invsigW + impM)

  wi <- genNormDataTweak(n, t(mw), Vw)
  # set factor variance to 1 to identify the model
  wi <- apply(wi, 2, function(x) x / sd(x))

  # sample phi for g-factor:
  phiw <- LaplacesDemon::rinvwishart(nu = n + p0w, S = t(wi) %*% (wi) + solve(R0w))

  return(list(psi = pp, lambda = ll, psiw = psiw, beta = beta, wi = wi, phiw = phiw))
}

drawStartMulti <- function(n, k, ns, pars, imat) {
  # measurement parameters
  invpsi <- rgamma(k, pars$a0k, pars$b0k)
  psi <- 1 / invpsi
  # structural parameters
  invpsiw <- rgamma(ns, pars$a0kw, pars$b0kw)
  psiw <- 1 / invpsiw

  # ------- factor scores for all factors ---------
  phiw <- LaplacesDemon::rinvwishart(nu = pars$p0w, S = (pars$R0w))

  wi <- MASS::mvrnorm(n, numeric(ns + 1), phiw)
  wi <- apply(wi, 2, function(x) x / sd(x))

  return(list(wi = wi, phiw = phiw))
}
