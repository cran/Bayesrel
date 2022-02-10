


omegaMultiB <- function(data, ns, n.iter, n.burnin, n.chains, thin, model, pairwise,
                        a0, b0, l0, A0, c0, d0, beta0, B0, p0, R0, param.out, callback, pbtick) {

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

  if (is.matrix(l0) || is.data.frame(l0)) {
    l0mat <- l0
  } else {
    l0mat <- matrix(0, k, ns)
    l0mat[imat] <- l0
  }

  beta0vec <- numeric(ns)
  beta0vec[1:ns] <- beta0

  pars <- list(H0k = rep(A0, ns), a0k = a0, b0k = b0, l0k = l0mat,
               H0kw = B0, a0kw = c0, b0kw = d0, beta0k = beta0vec,
               R0winv = diag(rep(R0, ns + 1)), p0w = p0)

  omsh <- matrix(0, n.chains, n.iter)
  omst <- matrix(0, n.chains, n.iter)
  impl_covs <- array(0, c(n.chains, n.iter, k, k))

  if (param.out) {
    lambdas <- array(0, c(n.chains, n.iter, k, ns))
    betas <- array(0, c(n.chains, n.iter, ns))
    thetas <- array(0, c(n.chains, n.iter, k))
    psis <- array(0, c(n.chains, n.iter, ns))
  }

  for (ai in 1:n.chains) {

    phiw <- diag(1 / rgamma(ns + 1, shape = pars$p0w / 2, scale = 2 / diag(pars$R0winv)))
    wi <- genNormDataLegit(n, numeric(ns + 1), phiw)

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
        oms <- omegasSeco(Lm, Bm, diag(params$theta), diag(c(1, params$psiw)))

        omsh[ai, i] <- oms[1]
        omst[ai, i] <- oms[2]

        cc <- implCovMulti(Lm, Bm, theta = diag(params$theta), psi = diag(c(1, params$psiw)))
        impl_covs[ai, i, , ] <- cc

        if (param.out) {
          lambdas[ai, i, , ] <- params$lambda
          betas[ai, i, ] <- params$beta
          thetas[ai, i, ] <- params$theta
          psis[ai, i, ] <- params$psiw
        }

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

        pbtick()
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
        oms <- omegasSeco(Lm, Bm, diag(params$theta), diag(c(1, params$psiw)))

        omsh[ai, i] <- oms[1]
        omst[ai, i] <- oms[2]

        impl_covs[ai, i, , ] <- implCovMulti(Lm, Bm, theta = diag(params$theta), psi = diag(c(1, params$psiw)))

        if (param.out) {
          lambdas[ai, i, , ] <- params$lambda
          betas[ai, i, ] <- params$beta
          thetas[ai, i, ] <- params$theta
          psis[ai, i, ] <- params$psiw
        }

        pbtick()
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

  if (param.out) {
    lambda_burn <- lambdas[, (n.burnin + 1):n.iter, , , drop = F]
    beta_burn <- betas[, (n.burnin + 1):n.iter, , drop = F]
    theta_burn <- thetas[, (n.burnin + 1):n.iter, , drop = F]
    psi_burn <- psis[, (n.burnin + 1):n.iter, , drop = F]

    lambda_out <- lambda_burn[, seq(1, dim(omt_burn)[2], thin), , , drop = F]
    beta_out <- beta_burn[, seq(1, dim(omt_burn)[2], thin), , drop = F]
    theta_out <- theta_burn[, seq(1, dim(omt_burn)[2], thin), , drop = F]
    psi_out <- psi_burn[, seq(1, dim(omt_burn)[2], thin), , drop = F]

    return(list(omh = omh_out, omt = omt_out, impl_covs = impl_covs_out, imputed_values = imputed,
                modfile = mod_opts, lambda = lambda_out, beta = beta_out, theta = theta_out, psi = psi_out))
  }

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

  R0winv <- pars$R0winv
  p0w <- pars$p0w

  ll <- matrix(0, k, ns)
  pp <- numeric(k)

  for (ii in 1:ns) {

    ids <- idex[[ii]]
    Ak <- solve(1 / H0k[ii] + t(wi[, ii + 1]) %*% wi[, ii + 1])
    ak <- Ak %*% (c(1 / H0k[ii]) %*% t(l0k[ids, ii]) + wi[, ii + 1] %*% data[, ids])
    bekk <- b0k + 0.5 * (t(data[, ids]) %*% data[, ids]
                         - t(ak) %*% solve(Ak) %*% ak
                         + (l0k[ids, ii] * (1 / H0k[ii])) %*% t(l0k[ids, ii]))
    bek <- diag(bekk)
    invpsi <- rgamma(length(ids), n / 2 + a0k, bek)
    psi <- 1 / invpsi
    lambda <- rnorm(length(ids), ak, sqrt(psi * as.vector(Ak)))

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
  beta <- rnorm(ns, akw * sqrt(phiw[1, 1]), sqrt(psiw * Akw))

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
  phiw <-  diag(1 / rgamma(ns + 1, shape = (n + pars$p0w) / 2, scale = 2 / diag(t(wi) %*% (wi) + pars$R0winv)))

  return(list(theta = pp, lambda = ll, psiw = psiw, beta = beta, wi = wi, phiw = phiw))
}

