
#'@export
print.strel <- function(x, ...) {
  if (!is.null(x$Bayes) & !is.null(x$freq)) {
    outB <- cbind(as.numeric(sprintf("%.5f", as.matrix(x$Bayes$est))),
                  as.numeric(sprintf("%.5f", as.matrix(x$Bayes$cred$low))),
                  as.numeric(sprintf("%.5f", as.matrix(x$Bayes$cred$up))))
    outF <- cbind(as.numeric(sprintf("%.5f", as.matrix(x$freq$est))),
                  as.numeric(sprintf("%.5f", as.matrix(x$freq$conf$low))),
                  as.numeric(sprintf("%.5f", as.matrix(x$freq$conf$up))))
    out <- rbind(outB, outF)
    rownames(out) <- c(names(x$Bayes$est), names(x$freq$est))
  }
  else if (!is.null(x$Bayes)) {
    out <- cbind(as.numeric(sprintf("%.5f", as.matrix(x$Bayes$est))),
                 as.numeric(sprintf("%.5f", as.matrix(x$Bayes$cred$low))),
                 as.numeric(sprintf("%.5f", as.matrix(x$Bayes$cred$up))))
    rownames(out) <- names(x$Bayes$est)
  }
  else if (!is.null(x$freq)) {
    out <- cbind(as.numeric(sprintf("%.5f", as.matrix(x$freq$est))),
                 as.numeric(sprintf("%.5f", as.matrix(x$freq$conf$low))),
                 as.numeric(sprintf("%.5f", as.matrix(x$freq$conf$up))))
    rownames(out) <- names(x$freq$est)

  } else {
    return(warning("no estimates calculated"))
  }
  colnames(out) <- c("point est", paste0(" ", x$interval * 100, " % CI lower"),
                     paste0(" ", x$interval * 100, "% CI upper"))
  cat("Call: \n")
  print.default(x$call)
  cat("\n")
  cat("Estimates of Single-Test Reliability Measures: \n")
  cat("\n")
  print.default(out)

}

#'@export
summary.strel <- function(object, ...) {

  out <- object

  class(out) <- "summary.strel"
  out
}

#'@export
print.summary.strel <- function(x, ...) {

  tmp_matrix <- list()
  if (!is.null(x$freq) && !is.null(x$Bayes)) {
    tmp_matrix$est <- rbind(as.matrix(x$Bayes$est),
                            as.matrix(x$freq$est))
    tmp_matrix$low <- rbind(as.matrix(x$Bayes$cred$low),
                                as.matrix(x$freq$conf$low))
    tmp_matrix$up <- rbind(as.matrix(x$Bayes$cred$up),
                               as.matrix(x$freq$conf$up))
  } else if (!is.null(x$Bayes)) {
    tmp_matrix$est <- as.matrix(x$Bayes$est)
    tmp_matrix$low <- as.matrix(x$Bayes$cred$low)
    tmp_matrix$up <- as.matrix(x$Bayes$cred$up)
  } else if (!is.null(x$freq)) {
    tmp_matrix$est <- as.matrix(x$freq$est)
    tmp_matrix$low <- as.matrix(x$freq$conf$low)
    tmp_matrix$up <- as.matrix(x$freq$conf$up)
  }

  n_row <- nrow(tmp_matrix$est)
  mat <- matrix(0, n_row, 3)
  mat[, 1] <- as.numeric(sprintf("%.5f", tmp_matrix$est))
  mat[, 2] <- as.numeric(sprintf("%.5f", tmp_matrix$low))
  mat[, 3] <- as.numeric(sprintf("%.5f", tmp_matrix$up))
  rownames(mat) <- rownames(tmp_matrix$est)
  colnames(mat) <- c("point est", paste0(" ", x$interval * 100, "% CI lower"),
                     paste0(" ", x$interval * 100, "% CI upper"))

  cat("Call: \n")
  print.default(x$call)
  cat("\n")
  cat("Results: \n")
  print(mat, right = FALSE)
  cat("\n")
  if (!is.null(x$n.iter)) {
    cat("Bayesian point est is the posterior mean \n")
  }
  if (length(grep("freq", rownames(tmp_matrix$est))) > 0) {

    if (!is.null(x$freq$inv.mat)) {
      cat("the number of bootstrap samples reduced to ")
      cat(x$freq$inv.mat)
      cat(" because some bootstrapped matrices were singular\n")
    }
    if ("omega" %in% x$estimates) {
      if (!is.null(x$freq$omega.pfa) & !is.null(x$freq$omega.error)) {
        cat("frequentist omega method is pfa\n")
        cat("omega confidence interval is estimated with bootstrap because the cfa did not find a solution\n")
      }
    }
  }
  if (!is.null(x$complete)) {
    cat("Missing data handling: using listwise deletion the number of complete cases is\n")
    cat(x$complete)
  }
  if (!is.null(x$miss_pairwise)) {
    cat("Missing data handling: using pairwise complete cases\n")
  }

  if (!is.null(x$Bayes$ifitem$est)) {
    n_row <- length(unlist(x$Bayes$ifitem$est[1])) + 1
    n_col <- 3
    names <- NULL
    for(z in 1:(n_row - 1)){
      names[z] <- paste0("x", z)
    }
    row_names <- c("original", names)

    for (i in seq_len(length(x$estimates))) {
      mat_ifitem_bay <- matrix(0, n_row, n_col)
      rownames(mat_ifitem_bay) <- row_names

      mat_ifitem_bay[1, ] <- as.numeric(sprintf("%.5f", c(unlist(tmp_matrix$est)[i],
                                                          unlist(tmp_matrix$low)[i],
                                                          unlist(tmp_matrix$up)[i])))
      mat_ifitem_bay[2:n_row, ] <- cbind(as.numeric(sprintf("%.5f", unlist(x$Bayes$ifitem$est[i]))),
                                         matrix(as.numeric(sprintf("%.5f", unlist(x$Bayes$ifitem$cred[i]))), n_row - 1, 2))
      colnames(mat_ifitem_bay) <- c("point est", paste0(" ", x$interval * 100, "% CI lower"),
                                    paste0(" ", x$interval * 100, "% CI upper"))
      cat("\n")
      cat(paste0("Bayesian ", x$estimate[i], " if item dropped: \n"))
      print.default(mat_ifitem_bay)
    }
  }

  if (!is.null(x$freq$ifitem)) {
    n_row <- length(unlist(x$freq$ifitem[1])) + 1
    n_col <- length(x$estimates)
    names <- NULL
    for (z in 1:(n_row-1)) {
      names[z] <- paste0("x", z)
    }
    row_names <- c("original", names)
    mat_ifitem_freq <- matrix(0, n_row, n_col)
    mat_ifitem_freq[1, ] <- as.numeric(sprintf("%.5f", unlist(tmp_matrix$est)[grep("freq", rownames(tmp_matrix$est))]))
    for (i in 1:n_col) {
      mat_ifitem_freq[2:n_row, i] <- as.numeric(sprintf("%.5f", unlist(x$freq$ifitem[i])))
    }
    colnames(mat_ifitem_freq) <- x$estimates
    row.names(mat_ifitem_freq) <- c("original", names)

    cat("\n")
    cat("Frequentist point estimate if item dropped: \n")
    print.default(mat_ifitem_freq)

    if ("omega" %in% x$estimates) {
      if (!is.null(x$omega.item.error)) {
        cat("frequentist omega method for item.dropped statistic is pfa because the cfa did not find a solution\n")
      }
    }
  }
}


#'@export
print.bomegas <- function(x, ...) {
  # prepare output matrix
  out <- rbind(as.numeric(sprintf("%.5f", c(x$omega_t$mean, x$omega_t$cred))),
               as.numeric(sprintf("%.5f", c(x$omega_h$mean, x$omega_h$cred))))
  rownames(out) <- c("omega_t", "omega_h")
  colnames(out) <- c("posterior mean", paste0(x$interval * 100, "% CI lower"),
                      paste0(x$interval * 100, "% CI upper"))

  # output:
  cat("Call: \n")
  print.default(x$call)
  cat("\n")
  print.default(out)
  if (x$listwise) {
    cat("\nComplete cases: ")
    cat(x$complete_cases)
  }
}

#'@export
print.omegasCFA <- function(x, ...) {
  # prepare output matrix
  out <- rbind(as.numeric(sprintf("%.5f", c(x$omega_t$est, x$omega_t$conf))),
               as.numeric(sprintf("%.5f", c(x$omega_h$est, x$omega_h$conf))))
  rownames(out) <- c("omega_t", "omega_h")
  colnames(out) <- c("point est", paste0(x$interval * 100, "% CI lower"),
                      paste0(x$interval * 100, "% CI upper"))

  # output:
  cat("Call: \n")
  print.default(x$call)
  cat("\n")
  print.default(out)
  if (x$listwise) {
    cat("\nComplete cases: ")
    cat(x$complete_cases)
  }

  if (!is.null(x$model$fit.measures)) {
    fit_names <- c("chisq", "df", "pvalue", "cfi", "tli",
                   "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue",
                   "aic", "bic", "usrmr", "usrmr.ci.lower", "usrmr.ci.upper", "usrmr.closefit.pvalue")
    measures <- unname(x$model$fit.measures[fit_names[1:11]])
    measures <- c(measures, x$model$srmr.summary[fit_names[12:15], ])
    measures <- as.numeric(sprintf("%.5f", measures))
    names(measures) <- fit_names
    measures <- as.data.frame(measures)
    cat("\nFit measures:\n")
    print(measures)
  }
}

#'@export
print.secoFit <- function(x, ...) {
  out <- list(LR = x$LR, BSRMR = x$srmr_pointEst, BRMSEA = x$rmsea_pointEst,
              BRMSEA_90_CI = x$rmsea_ci, BRMSEA_p.05 = x$p_rmsea)
  print(out)
}
