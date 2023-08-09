


omegaMultiBayes <- function(data, ns, n.iter, n.burnin, n.chains, thin, model, impute,
                        prior.params, param.out, callback, pbtick,
                        model.type) {

  n <- nrow(data)
  k <- ncol(data)

  model_opts <- indexMatrix(model, k, ns, colnames(data))
  # ---- get the index matrix aka which items load on which factors ----
  idex <- model_opts$idex
  imat <- model_opts$imat

  # ---- check missingsness ----
  inds <- which(is.na(data), arr.ind = TRUE)
  imputed <- NULL
  if (length(inds) > 0) {
    imputed <- array(0, c(n.chains, n.iter, nrow(inds)))
    inds_split <- split(inds[, 1], inds[, 2])
    unique_cols <- unique(inds[, 2])
  }


  # ---- sampling start --------

  if (is.matrix(prior.params[["l0"]]) || is.data.frame(prior.params[["l0"]])) {
    l0mat <- prior.params[["l0"]]
  } else {
    l0mat <- matrix(0, k, ns)
    l0mat[imat] <- prior.params[["l0"]]
  }

  if(model.type != "correlated") {
    if (model.type == "second-order") {
      beta0vec <- numeric(ns)
      beta0vec[1:ns] <- prior.params[["beta0"]]
    } else {
      beta0vec <- numeric(k)
      beta0vec[1:k] <- prior.params[["beta0"]]
    }
    R0winv <- diag(rep.int(prior.params[["R0"]], ns + 1))
    omsh <- matrix(0, n.chains, n.iter)

  } else {
    R0winv <- prior.params[["R0"]]
    beta0vec <- NULL
  }

  pars <- list(H0k = rep.int(prior.params[["A0"]], ns), a0k = prior.params[["a0"]], b0k = prior.params[["b0"]],
               l0k = l0mat, H0kw = prior.params[["B0"]], a0kw = prior.params[["c0"]], b0kw = prior.params[["d0"]],
               beta0k = beta0vec, R0winv = R0winv, p0w = prior.params[["p0"]])


  omst <- matrix(0, n.chains, n.iter)
  impl_covs <- array(0, c(n.chains, n.iter, k, k))

  if (param.out) {
    if (model.type != "correlated") {
      lambdas <- array(0, c(n.chains, n.iter, k, ns))
      betas <- array(0, c(n.chains, n.iter, length(beta0vec)))
      thetas <- array(0, c(n.chains, n.iter, k))
      psis <- array(0, c(n.chains, n.iter, ns))
    } else {
      lambdas <- array(0, c(n.chains, n.iter, k, ns))
      thetas <- array(0, c(n.chains, n.iter, k))
      phis <- array(0, c(n.chains, n.iter, ns, ns))
    }
  }

  ticks <- 1 # for progressbar
# ----------- model second-order ----------- #
  if (model.type == "second-order") {

    for (ai in 1:n.chains) {
      phiw <- diag(1 / rgamma(ns + 1, shape = pars$p0w / 2, scale = 2 / diag(pars$R0winv)))
      wi <- genNormDataLegit(n, numeric(ns + 1), phiw)

      if (impute) { # missing data
        dat_filled <- data
        dat_filled[inds] <- colMeans(data, na.rm = TRUE)[inds[, 2]]
        ms <- rep.int(0, k)

        for (i in 1:n.iter) {
          params <- sampleSecoParams(dat_filled, pars, wi, phiw, ns, idex, imat)
          wi <- params$wi
          phiw <- params$phiw
          # compute omega
          Lm <- cbind(0, params$lambda)
          Bm <- matrix(0, ns + 1, ns + 1)
          Bm[2:(ns + 1), 1] <- params$beta
          oms <- omegasSeco(Lm, Bm, diag(params$theta), diag(c(1, params$psiw)))

          omsh[ai, i] <- oms[1]
          omst[ai, i] <- oms[2]

          cc <- implCovMultiSeco(Lm, Bm, theta = diag(params$theta), psi = diag(c(1, params$psiw)))
          impl_covs[ai, i, , ] <- cc

          if (param.out) {
            lambdas[ai, i, , ] <- params$lambda
            betas[ai, i, ] <- params$beta
            thetas[ai, i, ] <- params$theta
            psis[ai, i, ] <- params$psiw
          }

          # substitute missing values one by one, where each value is drawn conditional on the rest of the data
          # see https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
          for (ic in seq_along(unique_cols)) {
            ccc <- unique_cols[[ic]]

            # for (ccc in cols) {
            rows <- inds[which(inds[, 2] == ccc), 1]
            mu1 <- ms[ccc]
            mu2 <- ms[-ccc]
            cc11 <- cc[ccc, ccc]
            cc21 <- cc[-ccc, ccc]
            cc12 <- cc[ccc, -ccc]
            cc22 <- cc[-ccc, -ccc]

            cc22_chol <- chol(cc22)
            ccq <- cc11 - Xt_invChol_X_2(cc22_chol, cc21)

            cc12_cc22_inv <- forwardsolve(cc22_chol, backsolve(cc22_chol, cc12, transpose = TRUE), upper.tri = TRUE)
            # == cc12 %*% chol2inv(cc22_chol) ==
            rows <- inds_split[[ic]]
            for (r in rows) {
              muq <- mu1 + cc12_cc22_inv %*% (as.numeric(dat_filled[r, -ccc]) - mu2)
              dat_filled[r, ccc] <- rnorm(1, muq, sqrt(ccq))
            }
          }
          imputed[ai, i, ] <- dat_filled[inds]

          ticks <- ticks + 1
          setTxtProgressBar(pbtick, ticks)
        }

      } else {
        for (i in 1:n.iter) {

          params <- sampleSecoParams(data, pars, wi, phiw, ns, idex, imat)
          wi <- params$wi
          phiw <- params$phiw

          # compute omega
          Lm <- cbind(0, params$lambda)
          Bm <- matrix(0, ns + 1, ns + 1)
          Bm[2:(ns + 1), 1] <- params$beta
          oms <- omegasSeco(Lm, Bm, diag(params$theta), diag(c(1, params$psiw)))

          omsh[ai, i] <- oms[1]
          omst[ai, i] <- oms[2]

          impl_covs[ai, i, , ] <- implCovMultiSeco(Lm, Bm, theta = diag(params$theta), psi = diag(c(1, params$psiw)))

          if (param.out) {
            lambdas[ai, i, , ] <- params$lambda
            betas[ai, i, ] <- params$beta
            thetas[ai, i, ] <- params$theta
            psis[ai, i, ] <- params$psiw
          }
          ticks <- ticks + 1
          setTxtProgressBar(pbtick, ticks)
        }
      }
    }
# ---------- model bi factor -------- #
  } else if (model.type == "bi-factor") {

    if (any(rowSums(model_opts$imat) > 1)) {
      stop("Crossloadings cannot be specified with the bi-factor model.")
    }

    for (ai in 1:n.chains) {
      phiw <- diag(1 / rgamma(ns + 1, shape = pars$p0w / 2, scale = 2 / diag(pars$R0winv)))
      wi <- genNormDataLegit(n, numeric(ns + 1), phiw)

      if (impute) { # missing data
        dat_filled <- data
        dat_filled[inds] <- colMeans(data, na.rm = TRUE)[inds[, 2]]
        ms <- rep.int(0, k)

        for (i in 1:n.iter) {
          params <- sampleBifParams(dat_filled, pars, wi, phiw, ns, idex, imat)
          wi <- params$wi
          phiw <- params$phiw
          # compute omega
          oms <- omegasBif(params$lambda, params$beta, diag(params$theta))
          omsh[ai, i] <- oms[1]
          omst[ai, i] <- oms[2]

          cc <- implCovMultiBif(params$lambda, params$beta, theta = diag(params$theta), psi = diag(params$psiw))
          impl_covs[ai, i, , ] <- cc

          if (param.out) {
            lambdas[ai, i, , ] <- params$lambda
            betas[ai, i, ] <- params$beta
            thetas[ai, i, ] <- params$theta
            psis[ai, i, ] <- params$psiw
          }

          # substitute missing values one by one, where each value is drawn conditional on the rest of the data
          # see https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
          for (ic in seq_along(unique_cols)) {
            ccc <- unique_cols[[ic]]
            #
            # for (ccc in cols) {
            rows <- inds[which(inds[, 2] == ccc), 1]
            mu1 <- ms[ccc]
            mu2 <- ms[-ccc]
            cc11 <- cc[ccc, ccc]
            cc21 <- cc[-ccc, ccc]
            cc12 <- cc[ccc, -ccc]
            cc22 <- cc[-ccc, -ccc]

            cc22_chol <- chol(cc22)
            ccq <- cc11 - Xt_invChol_X_2(cc22_chol, cc21)

            cc12_cc22_inv <- forwardsolve(cc22_chol, backsolve(cc22_chol, cc12, transpose = TRUE), upper.tri = TRUE)
            # == cc12 %*% chol2inv(cc22_chol) ==
            rows <- inds_split[[ic]]
            for (r in rows) {
              muq <- mu1 + cc12_cc22_inv %*% (as.numeric(dat_filled[r, -ccc]) - mu2)
              dat_filled[r, ccc] <- rnorm(1, muq, sqrt(ccq))
            }
          }
          imputed[ai, i, ] <- dat_filled[inds]

          ticks <- ticks + 1
          setTxtProgressBar(pbtick, ticks)
        }

      } else {
        for (i in 1:n.iter) {

          params <- sampleBifParams(data, pars, wi, phiw, ns, idex, imat)
          wi <- params$wi
          phiw <- params$phiw

          # compute omega
          oms <- omegasBif(params$lambda, params$beta, diag(params$theta))
          omsh[ai, i] <- oms[1]
          omst[ai, i] <- oms[2]

          impl_covs[ai, i, , ] <- implCovMultiBif(params$lambda, params$beta, theta = diag(params$theta), psi = diag(params$psiw))

          if (param.out) {
            lambdas[ai, i, , ] <- params$lambda
            betas[ai, i, ] <- params$beta
            thetas[ai, i, ] <- params$theta
            psis[ai, i, ] <- params$psiw
          }
          ticks <- ticks + 1
          setTxtProgressBar(pbtick, ticks)
        }
      }
    }
# ----------- model correlated ------------
  } else if (model.type == "correlated") {

    for (ai in 1:n.chains) {

      # draw starting values
      phiw <- LaplacesDemon::rinvwishart(nu = pars$p0w, S = pars$R0winv)
      wi <- genNormDataLegit(n, rep(0, ns), phiw)

      if (impute) { # missing data
        dat_filled <- data
        dat_filled[inds] <- colMeans(data, na.rm = TRUE)[inds[, 2]]
        ms <- rep.int(0, k)

        for (i in 1:n.iter) {

          params <- sampleCorrParams(dat_filled, pars, wi, phiw, ns, idex, imat)
          wi <- params$wi
          phiw <- params$phiw
          # compute omega
          omt <- omegaCorr(params$lambda, params$phiw, params$theta)

          omst[ai, i] <- omt

          cc <- implCovCorr(params$lambda, params$phiw, diag(params$theta))
          # diag(cc) <- diag(cc) + 1e-8
          impl_covs[ai, i, , ] <- cc

          if (param.out) {
            lambdas[ai, i, , ] <- params$lambda
            thetas[ai, i, ] <- params$theta
            phis[ai, i, , ] <- params$phiw
          }

          # substitute missing values one by one, where each value is drawn conditional on the rest of the data
          # see https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
          for (ic in seq_along(unique_cols)) {
            ccc <- unique_cols[[ic]]

            # for (ccc in cols) {
            rows <- inds[which(inds[, 2] == ccc), 1]
            mu1 <- ms[ccc]
            mu2 <- ms[-ccc]
            cc11 <- cc[ccc, ccc]
            cc21 <- cc[-ccc, ccc]
            cc12 <- cc[ccc, -ccc]
            cc22 <- cc[-ccc, -ccc]

            cc22_chol <- chol(cc22)
            ccq <- cc11 - Xt_invChol_X_2(cc22_chol, cc21)

            cc12_cc22_inv <- forwardsolve(cc22_chol, backsolve(cc22_chol, cc12, transpose = TRUE), upper.tri = TRUE)
            # == cc12 %*% chol2inv(cc22_chol) ==
            rows <- inds_split[[ic]]
            for (r in rows) {
              muq <- mu1 + cc12_cc22_inv %*% (as.numeric(dat_filled[r, -ccc]) - mu2)
              dat_filled[r, ccc] <- rnorm(1, muq, sqrt(ccq))
            }
          }

          imputed[ai, i, ] <- dat_filled[inds]

          ticks <- ticks + 1
          setTxtProgressBar(pbtick, ticks)
        }

      } else {
        for (i in 1:n.iter) {

          params <- sampleCorrParams(data, pars, wi, phiw, ns, idex, imat)
          wi <- params$wi
          phiw <- params$phiw
          # compute omega
          omt <- omegaCorr(params$lambda, params$phiw, params$theta)

          omst[ai, i] <- omt

          impl_covs[ai, i, , ] <- implCovCorr(params$lambda, params$phiw, diag(params$theta))

          if (param.out) {
            lambdas[ai, i, , ] <- params$lambda
            thetas[ai, i, ] <- params$theta
            phis[ai, i, , ] <- params$phiw
          }

          ticks <- ticks + 1
          setTxtProgressBar(pbtick, ticks)
        }
      }
    }
  } else {
    stop("Invalid model type entered.")
  }

  lambda_out <- NULL
  beta_out <- NULL
  theta_out <- NULL
  psi_out <- NULL
  phi_out <- NULL
  omh_out <- NULL

  if (model.type != "correlated") {

    omh_burn <- omsh[, (n.burnin + 1):n.iter, drop = FALSE]
    omh_out <- omh_burn[, seq(1, dim(omh_burn)[2], thin), drop = FALSE]

    if (param.out) {
      lambda_burn <- lambdas[, (n.burnin + 1):n.iter, , , drop = FALSE]
      beta_burn <- betas[, (n.burnin + 1):n.iter, , drop = FALSE]
      theta_burn <- thetas[, (n.burnin + 1):n.iter, , drop = FALSE]
      psi_burn <- psis[, (n.burnin + 1):n.iter, , drop = FALSE]

      lambda_out <- lambda_burn[, seq(1, dim(lambda_burn)[2], thin), , , drop = FALSE]
      beta_out <- beta_burn[, seq(1, dim(beta_burn)[2], thin), , drop = FALSE]
      theta_out <- theta_burn[, seq(1, dim(theta_burn)[2], thin), , drop = FALSE]
      psi_out <- psi_burn[, seq(1, dim(psi_burn)[2], thin), , drop = FALSE]
    }

  } else {

    if (param.out) {
      lambda_burn <- lambdas[, (n.burnin + 1):n.iter, , , drop = FALSE]
      theta_burn <- thetas[, (n.burnin + 1):n.iter, , drop = FALSE]
      phi_burn <- phis[, (n.burnin + 1):n.iter, , , drop = FALSE]

      lambda_out <- lambda_burn[, seq(1, dim(lambda_burn)[2], thin), , , drop = FALSE]
      theta_out <- theta_burn[, seq(1, dim(theta_burn)[2], thin), , drop = FALSE]
      phi_out <- phi_burn[, seq(1, dim(phi_burn)[2], thin), , , drop = FALSE]
    }
  }

  # burnin
  omt_burn <- omst[, (n.burnin + 1):n.iter, drop = FALSE]
  impl_covs_burn <- impl_covs[, (n.burnin + 1):n.iter, , , drop = FALSE]

  # thinning
  omt_out <- omt_burn[, seq(1, dim(omt_burn)[2], thin), drop = FALSE]
  impl_covs_out <- impl_covs_burn[, seq(1, dim(impl_covs_burn)[2], thin), , , drop = FALSE]


  return(list(omh = omh_out, omt = omt_out, impl_covs = impl_covs_out, imputed_values = imputed, modfile = model_opts,
              lambda = lambda_out, beta = beta_out, theta = theta_out, psi = psi_out, phi = phi_out))

}


sampleSecoParams <- function(data, pars, wi, phiw, ns, idex, imat) {

  n <- nrow(data)
  k <- ncol(data)
  # -------- measurement equation ------------
  H0k <- pars$H0k # prior multiplier matrix for lambdas variance (and covariance)
  l0k <- pars$l0k # prior lambdas
  a0k <- pars$a0k # prior shape parameter for gamma function for psis
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

  # first the crossloadings, if any
  crsl <- which(rowSums(imat) > 1)
  for (ii in crsl) {
    idFacs <- which(imat[ii, ])
    Ak <- solve(1 / H0k[idFacs] + t(wi[, idFacs + 1]) %*% wi[, idFacs + 1])
    ak <- Ak %*% (1 / H0k[idFacs] * (l0k[ii, idFacs]) + t(wi[, idFacs + 1]) %*% data[, ii])
    bekk <- b0k + 0.5 * (t(data[, ii]) %*% data[, ii]
                         - t(ak) %*% solve(Ak) %*% ak
                         + (t(l0k[ii, idFacs]) * (1 / H0k[idFacs])) %*% (l0k[ii, idFacs]))
    bek <- diag(bekk)
    invpsi <- rgamma(1, n / 2 + a0k, bek)
    psi <- 1 / invpsi
    lambda <- rnorm(length(idFacs), ak, sqrt(psi * diag(Ak)))
    ll[ii, idFacs] <- lambda
    pp[ii] <- psi
  }

  # the non-crossloadings
  for (iii in 1:ns) {
    idIts <- idex[[iii]][!(idex[[iii]] %in% crsl)]
    Ak_inv <- 1 / H0k[iii] + sum(wi[, iii + 1]^2)
    ak <- (c(1 / H0k[iii]) %*% t(l0k[idIts, iii]) + wi[, iii + 1] %*% data[, idIts]) / c(Ak_inv)
    # computes the diagonal of bekk directly - maybe precompute diag_Xt_X(data[, idIts] for all idIts?
    bekk <- b0k + 0.5 * (diag_Xt_X(data[, idIts]) - diag_Xt_X(ak) * Ak_inv + diag_X_Xt(l0k[idIts, iii, drop = FALSE]) / H0k[iii])
    bek <- bekk

    invpsi <- rgamma(length(idIts), n / 2 + a0k, bek)
    psi <- 1 / invpsi
    lambda <- rnorm(length(idIts), ak, sqrt(psi * as.vector(1 / Ak_inv)))

    ll[idIts, iii] <- lambda
    pp[idIts] <- psi
  }

  # ------- structural equation -----
  Akw <- 1 / (1 / H0kw + c(crossprod(wi[, 1])))
  akw <- Akw * (1 / H0kw * beta0k + crossprod(wi[, 1], wi[, 2:(ns + 1)]))

  # computes the diagonal of bekk directly
  bekkw <- b0kw + 0.5 * (diag_Xt_X(wi[, 2:(ns + 1)]) - diag_Xt_X(akw) / Akw + diag_Xt_X(matrix(beta0k, 1)) / H0kw)
  bekw <- bekkw

  invpsiw <- rgamma(ns, n / 2 + a0kw, bekw)
  psiw <- 1 / invpsiw
  beta <- rnorm(ns, akw * sqrt(phiw[1, 1]), sqrt(psiw * Akw))

  # in Lee it says to replace the usual inv Phi matrix when sampling the factor scores with a function
  # of the g-factor loadings and their residuals
  betaMat <- matrix(0, ns + 1, ns + 1)
  betaMat[2:(ns + 1), 1] <- beta
  ident <- diag(ns + 1)
  identInv <- ident + betaMat # note: (ident + betaMat) %*% (ident - betaMat) == ident

  invsigW <- crossprod(sweep(ident - betaMat, 1L, c(1, sqrt(psiw)), `/`))

  lll <- cbind(0, ll)
  lll_temp <- sweep(lll, 1L, pp, `/`) # solve(diag(pp)) %*% lll
  impM <- crossprod(lll_temp, lll)

  # This solve can be optimized using blockwise inversion, assuming that (1) impM is always diagonal and (2) invsigW is diagonal except for the first row/ column
  Vw <- solve(invsigW + impM)
  mw <- Vw %*% tcrossprod(t(lll_temp), data)

  wi <- genNormDataTweak(n, t(mw), Vw)
  # set factor variance to 1 to identify the model
  wi <- apply(wi, 2, function(x) x / sd(x))

  # sample phi for g-factor:
  scale <- 2 / (.colSums(wi * wi, nrow(wi), ncol(wi)) + diag(pars$R0winv))
  phiw <-  diag(1 / rgamma(ns + 1, shape = (n + pars$p0w) / 2, scale = scale))

  return(list(theta = pp, lambda = ll, psiw = psiw, beta = beta, wi = wi, phiw = phiw))
}



sampleBifParams <- function(data, pars, wi, phiw, ns, idex, imat) {

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

  R0winv <- pars$R0winv
  p0w <- pars$p0w

  ll <- matrix(0, k, ns)
  pp <- numeric(k)

  H0 <- c(H0kw, H0k)
  beta0 <- cbind(beta0k, l0k)
  Aka_inv <- diag(1/(H0)) + crossprod(wi)
  Aka <- solve(Aka_inv)
  aka <- Aka %*% (1/H0 * t(beta0)  + crossprod(wi, data))
  bekka <- b0k + 0.5 * (crossprod(data) - crossprod(aka, Aka_inv) %*% aka + tcrossprod(beta0 * (1/H0), beta0))
  beka <- diag(bekka)

  invpsi_a <- rgamma(k, n / 2 + a0k, beka)
  invPsi <- diag(invpsi_a)
  psi_a <- 1/invpsi_a

  lll <- matrix(0, k, ns + 1)
  pp <- psi_a

  imatb <- cbind(TRUE, imat)
  mms <- (t(aka) %*% sqrt(phiw))[imatb]
  dAka <- diag(Aka)
  lll[imatb] <- rnorm(length(mms), mms,
                      sqrt(c(psi_a * dAka[1], psi_a * rep.int(dAka[2:(ns + 1)], times = colSums(imat)))))

  invM <- diag(ns + 1)
  tmpM <- t(lll) %*% invPsi %*% lll
  mw <- solve(invM + tmpM) %*% t(lll) %*% invPsi %*% t(data)
  Vw <- solve(invM + tmpM)

  beta <- lll[, 1]
  ll <- lll[, -1]
  psiw <- rep(1, ns + 1)

  wi <- genNormDataTweak(n, t(mw), Vw)
  # set factor variance to 1 to identify the model
  wi <- apply(wi, 2, function(x) x / sd(x))

  # sample phi for g-factor:
  scale <- 2 / (.colSums(wi * wi, nrow(wi), ncol(wi)) + diag(pars$R0winv))
  phiw <-  diag(1 / rgamma(ns + 1, shape = (n + pars$p0w) / 2, scale = scale))
  # phiw <-  diag(1 / rgamma(ns + 1, shape = (n + pars$p0w) / 2, scale = 2 / diag(t(wi) %*% (wi) + pars$R0winv)))

  return(list(theta = pp, lambda = ll, psiw = psiw, beta = beta, wi = wi, phiw = phiw))
}


sampleCorrParams <- function(data, pars, wi, phiw, ns, idex, imat) {

  n <- nrow(data)
  k <- ncol(data)

  H0k <- pars$H0k # prior multiplier for lambdas variance
  l0k <- pars$l0k # prior lambdas
  a0k <- pars$a0k # prior shape parameter for gamma function for psis
  b0k <- pars$b0k # prior rate parameter for gamma for psi

  p0 <- pars$p0w
  R0inv <- pars$R0winv

  invpsi <- rgamma(k, a0k, b0k)
  psi <- 1/invpsi
  invPsi <- diag(invpsi)

  ll <- matrix(0, k, ns)
  pp <- numeric(k)

# start with the crossloadings
  crsl <- which(rowSums(imat) > 1)
  for (ii in crsl) {
    idFacs <- which(imat[ii, ])
    Ak <- solve(1 / H0k[idFacs] + t(wi[, idFacs]) %*% wi[, idFacs])
    ak <- Ak %*% (1 / H0k[idFacs] * (l0k[ii, idFacs]) + t(wi[, idFacs]) %*% data[, ii])
    bekk <- b0k + 0.5 * (t(data[, ii]) %*% data[, ii]
                         - t(ak) %*% solve(Ak) %*% ak
                         + (t(l0k[ii, idFacs]) * (1 / H0k[idFacs])) %*% (l0k[ii, idFacs]))
    bek <- diag(bekk)
    invpsi <- rgamma(1, n / 2 + a0k, bek)
    psi <- 1 / invpsi
    lambda <- rnorm(length(idFacs), ak, sqrt(psi * diag(Ak)))

    ll[ii, idFacs] <- lambda
    pp[ii] <- psi
  }

# then the non crossloadings
  for (iii in 1:ns) {
    idIts <- idex[[iii]][!(idex[[iii]] %in% crsl)]
    Ak_inv <- 1 / H0k[iii] + sum(wi[, iii]^2)
    ak <- (c(1 / H0k[iii]) %*% t(l0k[idIts, iii]) + wi[, iii] %*% data[, idIts]) / c(Ak_inv)
    # computes the diagonal of bekk directly - maybe precompute diag_Xt_X(data[, idIts] for all idIts?
    bekk <- b0k + 0.5 * (diag_Xt_X(data[, idIts]) - diag_Xt_X(ak) * Ak_inv + diag_X_Xt(l0k[idIts, iii, drop = FALSE]) / H0k[iii])
    bek <- bekk

    invpsi <- rgamma(length(idIts), n / 2 + a0k, bek)
    psi <- 1 / invpsi
    lambda <- rnorm(length(idIts), ak, sqrt(psi * as.vector(1 / Ak_inv)))

    ll[idIts, iii] <- lambda
    pp[idIts] <- psi
  }

  phiw <- cov2cor(phiw) # ensures that the factor variances are 1
  invPhi <- solve(phiw)
  invPsi <- diag(1/pp)
  # sample wi posterior:
  crp <- crossprod(ll, invPsi)
  m <- solve(invPhi + crp %*% ll) %*% t(ll) %*% invPsi %*% t(data)
  V <- solve(invPhi + crp %*% ll)

  wi <- genNormDataTweak(n, t(m), V)
  # set factor variance to 1 to identify the model
  # wi <- apply(wi, 2, function(x) x/sd(x))

  # sample phi:
  phiw <- LaplacesDemon::rinvwishart(nu = n + p0, S = t(wi) %*% (wi) + R0inv)

  return(list(lambda = ll, theta = pp, phiw = phiw, wi = wi))

}

