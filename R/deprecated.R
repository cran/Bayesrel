

#################################################
#------------- Deprecated functions -------------#
#################################################


#' @title Deprecated functions in package \pkg{Bayesrel}.
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality are also mentioned. Help pages for deprecated functions are
#'   available at \code{help("<function>-deprecated")}.
#' @name Bayesrel-deprecated
#' @keywords internal
NULL


#' graphical posterior predictive check for the 1-factor omega model,
#' based on covariance matrix eigenvalues
#'
#' @description
#' gives posterior predictive check for the 1-factor model:
#' comparison between model implied covariance matrix and sample covariance matrix
#' also displays frequentist fit indices
#'
#' @param x A strel output object (list)
#' @param data A matrix or data.frame containing the data set that produced x
#' @param ppc A logical indicating if the PPC should be printed or not,
#' the default is TRUE
#' @param cutoff A value to compare the posterior sample of RMSEAs against.
#' The result will contain the probability that the RMSEA is smaller than the
#' cutoff value
#' @param ci A value between 0 and 1 indicating the credible interval for the RMSEA
#'
#' @references{
#' \insertRef{Garnier-Villarreal2020AdaptingFitIndices}{Bayesrel}
#' }
#'
#' @importFrom stats complete.cases
#'
#' @name omega_fit-deprecated
#' @usage omega_fit(x, data, ppc = TRUE, cutoff = .08, ci = .90)
#' @seealso \code{\link{Bayesrel-deprecated}}
#' @keywords internal
NULL
#'
#' @rdname Bayesrel-deprecated
#' @section \code{omega_fit}:
#' For \code{omega_fit}, use \code{\link{omegaFit}}.
#'
#' @export
omega_fit <- function(x, data, ppc = TRUE, cutoff = .08, ci = .90) {
  .Deprecated("omegaFit")

  if (!("omega" %in% x$estimates)) {
    return(warning("please run the analysis with omega as an estimate"))
  }
  if (!is.null(x$freq$omega.pfa)) {
    warning("cannot compute fit.measures for the pfa method")
  }

  data <- scale(data, scale = FALSE)

  outF <- outB <- NULL

  if (!is.null(x$Bayes)) {
    if (!is.null(x$miss_pairwise)) {
      sigma <- cov(data, use = "pairwise.complete.obs")
    } else {
      sigma <- cov(data, use = "complete.obs")
    }

    lambda <- x$Bayes$loadings
    psi <- x$Bayes$resid_var
    phi <- x$Bayes$f_var
    nsamp <- nrow(lambda)
    implieds <- array(NA, c(nsamp, ncol(lambda), ncol(lambda)))

    for (i in 1:nsamp) {
      implieds[i, , ] <- lambda[i, ] %*% t(phi[i]) %*% t(lambda[i, ]) + diag(psi[i, ])
    }

    if (ppc) {
      ee_impl <- matrix(0, nsamp, ncol(sigma))
      for (i in 1:nsamp) {
        dtmp <- MASS::mvrnorm(nrow(data), rep(0, ncol(sigma)), implieds[i, , ])
        ee_impl[i, ] <- eigen(cov(dtmp), only.values = TRUE)$values
      }
      qq_ee_low <- apply(ee_impl, 2, quantile, prob = .025)
      qq_ee_up <- apply(ee_impl, 2, quantile, prob = .975)

      ymax <- max(ee_impl, eigen(sigma, only.values = TRUE)$values)

      plot(eigen(sigma)$values, axes = F, ylim = c(0, ymax), ylab = "Eigenvalue", xlab = "Eigenvalue No.",
           pch = 8, cex = 0)

      arrows(x0 = seq_len(ncol(sigma)), x1 = seq_len(ncol(sigma)), y0 = qq_ee_low, y1 = qq_ee_up,
             col = "gray55", angle = 90, code = 3, length = .06, lwd = 2.5)

      lines(eigen(sigma)$values, type = "p", pch = 8, cex =.7)
      axis(side = 1, at = seq_len(ncol(sigma)))
      axis(side = 2, las = 1)

      legend(ncol(sigma) / 3 * 1.1, ymax * (2 / 3),
             legend = c("Dataset Covariance Matrix", "Model Implied Covariance Matrices"),
             col=c("black", "gray50"), box.lwd = .7, lty = c(0, 1), lwd = c(1, 2.5), pch = c(8, 30), cex = .8)
    }

    #### fit indices ###
    n <- nrow(data)
    k <- ncol(data)
    pstar <- k * (k + 1) / 2 # unique elements in the covariance matrix, variances + covariances

    ### Chisqs ###
    LL1 <- sum(dmultinorm(data, sigma)) # loglikelihood saturated model
    LR_obs <- apply(implieds, 1, LRblav, data = data, basell = LL1) # loglikelihoods tested model

    lsm <- colMeans(lambda)
    tsm <- colMeans(psi)
    pm <- mean(phi)
    implM <- lsm %*% t(1) %*% t(lsm) + diag(tsm) # mean model implied matrix, factor variance set to 1
    # implM <- apply(implieds, c(2, 3), mean)
    Dtm <- LRblav(data, implM, LL1) # deviance of the mean model implied cov matrix

    Dm <- mean(LR_obs) # mean deviance of the model implied cov matrices
    pD <- Dm - Dtm # effective number of parameters (free parameters)
    rmsea <- BRMSEA(LR_obs, pstar, pD, n)
    srmr_m <- SRMR(sigma, implM)
    # srmr <- apply(implieds, 1, SRMR, cdat = sigma)

    prob <- mean(rmsea < cutoff)

    rmsea_ci <- as.numeric(coda::HPDinterval(coda::mcmc(rmsea), prob = ci))
    names(rmsea_ci) <- c("lower", "upper")

    outB <- list(LR = Dtm, srmr = srmr_m, rmsea = mean(rmsea),
                 rmsea_ci = rmsea_ci, p_rmsea = prob)
  }

  if (!is.null(x$freq$omega_fit)) {
    outF <- list(freq_fit = x$freq$omega_fit)
  }

  return(list(Bayes = outB, freq = outF))
}


#' model fit for the second-order factor model,
#' @description
#' Fit indices and posterior predictive check for the higher-factor model:
#' comparison between posterior sample of model implied covariance matrices
#' and sample covariance matrix. Gray bars should enclose the black dots for good fit.
#' Also prints fit indices, LR (likelihood-ratio), RMSEA, SRMR.
#' The RMSEA is from Garnier-Villareal & Jorgensen (2020)
#'
#' @param x A bomegas output object (list)
#' @param data A matrix or data.frame containing the data set that produced x
#' @param ppc A logical indicating if the PPC should be printed or not,
#' the default is TRUE
#' @param cutoff A value to compare the posterior sample of RMSEAs against.
#' The result will contain the probability that the RMSEA is smaller than the
#' cutoff value
#' @param ci A value between 0 and 1 indicating the credible interval for the RMSEA
#'
#'
#' @references{
#' \insertRef{Garnier-Villarreal2020AdaptingFitIndices}{Bayesrel}
#' }
#'
#' @importFrom stats complete.cases
#'
#' @name seco_fit-deprecated
#' @usage seco_fit(x, data, ppc = TRUE, cutoff = .08, ci = .90)
#' @seealso \code{\link{Bayesrel-deprecated}}
#' @keywords internal
NULL
#'
#' @rdname Bayesrel-deprecated
#' @section \code{seco_fit}:
#' For \code{seco_fit}, use \code{\link{multiFit}}.
#'
#' @export
seco_fit <- function(x, data, ppc = TRUE, cutoff = .08, ci = .90) {
  .Deprecated("multiFit")

  data <- scale(data, scale = FALSE)

  if (x$pairwise) {
    sigma <- cov(data, use = "pairwise.complete.obs")
    n.data <- nrow(data)
  } else {
    sigma <- cov(data, use = "complete.obs")
    n.data <- nrow(data[complete.cases(data),])
  }

  if (ppc) {
    # PPC plot
    nsamp <- dim(x$implCovs)[1]
    ee_impl <- matrix(0, nsamp, ncol(sigma))
    for (i in 1:nsamp) {
      ctmp <-x$implCovs[i, , ]
      dtmp <- MASS::mvrnorm(n.data, rep(0, ncol(sigma)), ctmp)
      ee_impl[i, ] <- eigen(cov(dtmp), only.values = TRUE)$values
    }
    qq_ee_low <- apply(ee_impl, 2, quantile, prob = .025)
    qq_ee_up <- apply(ee_impl, 2, quantile, prob = .975)

    ymax <- max(ee_impl, eigen(sigma, only.values = TRUE)$values)

    par(mar=c(5.1, 4.5, 0.7, 2.1))
    plot(eigen(sigma)$values, axes = F, ylim = c(0, ymax), pch = 8, cex = 0, xlab = "", ylab = "")

    arrows(x0 = seq_len(ncol(sigma)), x1 = seq_len(ncol(sigma)), y0 = qq_ee_low, y1 = qq_ee_up,
           col = "gray55", angle = 90, code = 3, length = .06, lwd = 2.5)

    lines(eigen(sigma)$values, type = "p", pch = 8, cex =.7)
    axis(side = 1, at = seq_len(ncol(sigma)), cex.axis = 1.4)
    axis(side = 2, las = 1, cex.axis = 1.4)
    title(xlab = "Eigenvalue No.", ylab = "Eigenvalue", cex.lab = 1.4)

    legend(ncol(sigma) / 3 * 1.1, ymax  *(2 / 3),
           legend = c("Dataset Covariance Matrix", "Model Implied Covariance Matrices"),
           col=c("black", "gray50"), box.lwd = .7, lty = c(0, 1), lwd = c(1, 2.5), pch = c(8, 0), cex = 1.2,
           pt.cex = c(1, 0))
  }

  #### fit indices ###
  n <- n.data
  k <- ncol(data)
  pstar <- k * (k + 1) / 2 # unique elements in the covariance matrix, variances + covariances
  implieds <- x$implCovs

  ### Chisqs ###
  LL1 <- sum(dmultinorm(data, sigma)) # loglikelihood saturated model
  LR_obs <- apply(implieds, 1, LRblav, data = data, basell = LL1) # loglikelihoods tested model

  implM <- apply(implieds, c(2, 3), mean) # mean model implied matrix
  Dtm <- LRblav(data, implM, LL1) # deviance of the mean model implied cov matrix

  Dm <- mean(LR_obs) # mean deviance of the model implied cov matrices
  pD <- Dm - Dtm # effective number of parameters (free parameters)
  rmsea <- BRMSEA(LR_obs, pstar, pD, n)
  srmr_m <- SRMR(sigma, implM)
  # srmr <- apply(implieds, 1, SRMR, cdat = sigma)

  prob <- mean(rmsea < cutoff)

  rmsea_ci <- as.numeric(coda::HPDinterval(coda::mcmc(rmsea), prob = ci))
  names(rmsea_ci) <- paste0(ci*100, "% ",  c("lower", "upper"))
  out <- list(LR = Dtm, srmr = srmr_m, rmsea = mean(rmsea),
              rmsea_ci = rmsea_ci, p_rmsea = prob)

  return(out)
}



#' prior and posterior probability of estimate being bigger than threshold
#' @description
#' takes a mcmc posterior sample of any of the single test reliability estimates
#' and calculates the prior and posterior probability of the estimate being bigger
#' or smaller than an arbitrary value (priors are stored in the package)
#'
#' @param x A strel output object (list)
#' @param estimate A character string indicating what estimate to plot from the strel output object
#' @param low.bound A number for the threshold to be tested against
#'
#' @name p_strel-deprecated
#' @usage p_strel(x, estimate, low.bound)
#' @seealso \code{\link{Bayesrel-deprecated}}
#' @keywords internal
NULL
#'
#' @rdname Bayesrel-deprecated
#' @section \code{p_strel}:
#' For \code{p_strel}, use \code{\link{pStrel}}.
#'
#' @export
p_strel <- function(x, estimate, low.bound) {
  .Deprecated("pStrel")

  posi1 <- grep(estimate, x$estimates, ignore.case = TRUE)
  samp <- as.vector(x$Bayes$samp[[posi1]])
  obj <- ecdf(samp)
  post_prob <- 1 - obj(low.bound)

  # prior prob
  n.item <- dim(x$Bayes$covsamp)[3]

  if (is.null(x$priors$df0)) {
    x$priors$df0 <- n.item
  }

  prior <- priorSampUni(n.item, estimate, k0 = x$priors$k0, df0 = x$priors$df0, a0 = x$priors$a0, b0 = x$priors$b0,
                        m0 = x$priors$m0)

  end <- length(prior[["x"]])
  poslow <- end - sum(prior[["x"]] > low.bound)
  prior_prob <- sum(prior[["y"]][poslow:end]) / sum(prior[["y"]])
  out <- c(prior_prob, post_prob)
  names(out) <- c("prior_prob", "posterior_prob")
  return(out)
}



#' prior and posterior probability of omega_t and omega_h being bigger than thresholds
#' @description
#' takes mcmc posterior samples of omega_t and omega_h
#' and calculates the prior and posterior probability of the estimate being bigger
#' or smaller than an arbitrary value
#'
#' @param x A strel output object (list)
#' @param cutoff.t A number indicating the threshold for omega_t
#' @param cutoff.h A number indicating the threshold for omega_h
#'
#' @name p_omegas-deprecated
#' @usage p_omegas(x, cutoff.t = .80, cutoff.h = .60)
#' @seealso \code{\link{Bayesrel-deprecated}}
#' @keywords internal
NULL
#'
#' @rdname Bayesrel-deprecated
#' @section \code{p_omegas}:
#' For \code{p_omegas}, use \code{\link{pOmegas}}.
#'
#' @export
p_omegas <- function(x, cutoff.t = .80, cutoff.h = .60) {
  .Deprecated("pOmegas")

  sampt <- as.vector(x$omega_t$chains)
  samph <- as.vector(x$omega_h$chains)
  objt <- ecdf(sampt)
  objh <- ecdf(samph)
  post_prob_t <- 1 - objt(cutoff.t)
  post_prob_h <- 1 - objh(cutoff.h)

  # prior prob
  priors <- omegasSecoPrior(x, nsamp = 2e3)

  priort <- ecdf(priors$omt_prior)
  priorh <- ecdf(priors$omh_prior)

  prior_prob_t <- 1 - priort(cutoff.t)
  prior_prob_h <- 1 - priorh(cutoff.h)

  out <- matrix(c(prior_prob_t, prior_prob_h, post_prob_t, post_prob_h, cutoff.t, cutoff.h),
                3, 2, byrow = TRUE)
  colnames(out) <- c("omega_t", "omega_h")
  rownames(out) <- c("prior_prob", "posterior_prob", "cutoff")
  return(out)
}


#' model fit for the second-order factor model,
#' @description
#' Fit indices and posterior predictive check for the higher-factor model:
#' comparison between posterior sample of model implied covariance matrices
#' and sample covariance matrix. Gray bars should enclose the black dots for good fit.
#' Also prints fit indices, LR (likelihood-ratio), RMSEA, SRMR.
#' The RMSEA is from Garnier-Villareal & Jorgensen (2020)
#'
#' @param x A bomegas output object (list)
#' @param data A matrix or data.frame containing the data set that produced x
#' @param ppc A logical indicating if the PPC should be printed or not,
#' the default is TRUE
#' @param cutoff A value to compare the posterior sample of RMSEAs against.
#' The result will contain the probability that the RMSEA is smaller than the
#' cutoff value
#' @param ci A value between 0 and 1 indicating the credible interval for the RMSEA
#'
#'
#' @references{
#' \insertRef{Garnier-Villarreal2020AdaptingFitIndices}{Bayesrel}
#' }
#'
#' @importFrom stats complete.cases
#'
#' @name secoFit-deprecated
#' @usage secoFit(x, data, ppc = TRUE, cutoff = .08, ci = .90)
#' @seealso \code{\link{Bayesrel-deprecated}}
#' @keywords internal
NULL
#'
#' @rdname Bayesrel-deprecated
#' @section \code{secoFit}:
#' For \code{secoFit}, use \code{\link{multiFit}}.
#' @export
secoFit <- function(x, data, ppc = TRUE, cutoff = .08, ci = .90) {
  .Deprecated("multiFit")

  data <- scale(data, scale = FALSE)

  if (x$pairwise) {
    sigma <- cov(data, use = "pairwise.complete.obs")
    n.data <- nrow(data)
  } else {
    sigma <- cov(data, use = "complete.obs")
    n.data <- nrow(data[complete.cases(data),])
  }

  if (ppc) {
    # PPC plot
    nsamp <- dim(x$implCovs)[1]
    ee_impl <- matrix(0, nsamp, ncol(sigma))
    for (i in 1:nsamp) {
      ctmp <-x$implCovs[i, , ]
      dtmp <- MASS::mvrnorm(n.data, rep(0, ncol(sigma)), ctmp)
      ee_impl[i, ] <- eigen(cov(dtmp), only.values = TRUE)$values
    }
    qq_ee_low <- apply(ee_impl, 2, quantile, prob = .025)
    qq_ee_up <- apply(ee_impl, 2, quantile, prob = .975)

    ymax <- max(ee_impl, eigen(sigma, only.values = TRUE)$values)

    par(mar=c(5.1, 4.5, 0.7, 2.1))
    plot(eigen(sigma)$values, axes = F, ylim = c(0, ymax), pch = 8, cex = 0, xlab = "", ylab = "")

    arrows(x0 = seq_len(ncol(sigma)), x1 = seq_len(ncol(sigma)), y0 = qq_ee_low, y1 = qq_ee_up,
           col = "gray55", angle = 90, code = 3, length = .06, lwd = 2.5)

    lines(eigen(sigma)$values, type = "p", pch = 8, cex =.7)
    axis(side = 1, at = seq_len(ncol(sigma)), cex.axis = 1.4)
    axis(side = 2, las = 1, cex.axis = 1.4)
    title(xlab = "Eigenvalue No.", ylab = "Eigenvalue", cex.lab = 1.4)

    legend(ncol(sigma) / 3 * 1.1, ymax  *(2 / 3),
           legend = c("Dataset Covariance Matrix", "Model Implied Covariance Matrices"),
           col=c("black", "gray50"), box.lwd = .7, lty = c(0, 1), lwd = c(1, 2.5), pch = c(8, 0), cex = 1.2,
           pt.cex = c(1, 0))
  }

  #### fit indices ###
  n <- n.data
  k <- ncol(data)
  pstar <- k * (k + 1) / 2 # unique elements in the covariance matrix, variances + covariances
  implieds <- x$implCovs

  ### Chisqs ###
  LL1 <- sum(dmultinorm(data, sigma)) # loglikelihood saturated model
  LR_obs <- apply(implieds, 1, LRblav, data = data, basell = LL1) # loglikelihoods tested model

  implM <- apply(implieds, c(2, 3), mean) # mean model implied matrix
  Dtm <- LRblav(data, implM, LL1) # deviance of the mean model implied cov matrix

  Dm <- mean(LR_obs) # mean deviance of the model implied cov matrices
  pD <- Dm - Dtm # effective number of parameters (free parameters)
  rmsea <- BRMSEA(LR_obs, pstar, pD, n)
  srmr_m <- SRMR(sigma, implM)
  srmr <- apply(implieds, 1, SRMR, cdat = sigma)

  prob <- mean(rmsea < cutoff)

  rmsea_ci <- as.numeric(coda::HPDinterval(coda::mcmc(rmsea), prob = ci))
  names(rmsea_ci) <- paste0(ci*100, "% ",  c("lower", "upper"))
  out <- list(LR = Dtm, srmr_pointEst = srmr_m, srmr_samp = srmr, rmsea_pointEst = mean(rmsea),
              rmsea_ci = rmsea_ci, p_rmsea = prob, rmsea_samp = rmsea)

  class(out) <- "secoFit"
  return(out)
}


