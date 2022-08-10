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
#' @examples omegaFit(strel(asrm, "omega", n.chains = 2, n.iter = 100), data = asrm)
#'
#' @references{
#' \insertRef{Garnier-Villarreal2020AdaptingFitIndices}{Bayesrel}
#' }
#' @importFrom stats complete.cases
#'
#' @export
omegaFit <- function(x, data, ppc = TRUE, cutoff = .08, ci = .90) {

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
#' @examples secoFit(bomegas(upps, n.factors = 5, n.chains = 2, n.iter = 100,
#' n.burnin = 50, missing = "listwise"), upps)
#'
#' @references{
#' \insertRef{Garnier-Villarreal2020AdaptingFitIndices}{Bayesrel}
#' }
#'
#' @importFrom stats complete.cases
#' @export
secoFit <- function(x, data, ppc = TRUE, cutoff = .08, ci = .90) {

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


