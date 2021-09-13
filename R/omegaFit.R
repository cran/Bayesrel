
#' graphical posterior predictive check for the 1-factor omega model,
#' based on covariance matrix eigenvalues
#'
#' @description
#' gives posterior predictive check for the 1-factor model:
#' comparison between model implied covariance matrix and sample covariance matrix
#' also displays frequentist fit indices
#'
#' @param x A strel output object (list)
#'
#' @examples omega_fit(strel(asrm, "omega", n.chains = 2, n.iter = 100))
#'
#' @export
#'
omega_fit <- function(x) {
  if (!("omega" %in% x$estimates)) {
    return(warning("please run the analysis with omega as an estimate"))
  }
  if (!is.null(x$freq$omega.pfa)) {
    warning("cannot compute fit.measures for the pfa method")
  }

  if (!is.null(x$Bayes)) {
    if (!is.null(x$miss_pairwise)) {
      sigma <- cov(x$data, use = "pairwise.complete.obs")
    } else {
      sigma <- cov(x$data, use = "complete.obs")
    }

    lambda <- x$Bayes$loadings
    psi <- x$Bayes$resid_var

    nsamp <- nrow(lambda)
    ee_impl <- matrix(0, nsamp, ncol(sigma))
    for (i in 1:nsamp) {
      ctmp <- lambda[i, ] %*% t(lambda[i, ]) + diag(psi[i, ])
      dtmp <- MASS::mvrnorm(nrow(x$data), rep(0, ncol(sigma)), ctmp)
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
  if (!is.null(x$freq$omega_fit)){
    return(x$freq$omega_fit)}

}


#' graphical posterior predictive check for the higher-order factor model,
#'
#' @description
#' gives posterior predictive check for the higher-factor model:
#' comparison between posterior sample of model implied covariance matrices
#' and sample covariance matrix. Gray bars should enclose the black dots for good fit.
#'
#' @param x A bomegas output object (list)
#' @param data A matrix or data.frame containing the data set that produced x
#'
#' @examples seco_fit(bomegas(upps, n.factors = 5, n.chains = 2, n.iter = 100,
#' n.burnin = 50, missing = "listwise"), upps)
#'
#' @importFrom stats complete.cases
#' @export
seco_fit <- function(x, data) {

  if (x$pairwise) {
    sigma <- cov(data, use = "pairwise.complete.obs")
    n.data <- nrow(data)
  } else {
    sigma <- cov(data, use = "complete.obs")
    n.data <- nrow(data[complete.cases(data),])
  }

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
