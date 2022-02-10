
#'@export
print.strel <- function(x, ...) {
  if (!is.null(x$Bayes) & !is.null(x$freq)) {
    out <- cbind(as.data.frame(as.matrix(x$Bayes$cred$low)), as.data.frame(as.matrix(x$Bayes$cred$up)),
                 as.data.frame(as.matrix(x$freq$conf$low)), as.data.frame(as.matrix(x$freq$conf$up)))
    colnames(out) <- c("Bay.lower", "Bay.upper","freq.lower", "freq.upper")
  }
  else if (!is.null(x$Bayes)) {
    out <- cbind(as.data.frame(as.matrix(x$Bayes$cred$low)), as.data.frame(as.matrix(x$Bayes$cred$up)))
    colnames(out) <- c("Bay.lower", "Bay.upper")
  }
  else if (!is.null(x$freq)) {
    out <- cbind(as.data.frame(as.matrix(x$freq$conf$low)), as.data.frame(as.matrix(x$freq$conf$up)))
    colnames(out) <- c("freq.lower", "freq.upper")
  } else {
    return(warning("no estimates calculated"))
  }
  row.names(out) <- x$estimates
  cat("Call: \n")
  print.default(x$call)
  cat("\n")
  cat(x$interval,"% Interval Estimates of Single-Test Reliability Measures: \n")
  cat("\n")
  print(out)

}

#'@export
summary.strel <- function(object, ...) {

  out_matrix <- list()
  if (!is.null(object$freq) & !is.null(object$Bayes)) {
    out_matrix$est <- rbind(as.data.frame(as.matrix(object$Bayes$est)),
                            as.data.frame(as.matrix(object$freq$est)))
    out_matrix$int$low <- rbind(as.data.frame(as.matrix(object$Bayes$cred$low)),
                                as.data.frame(as.matrix(object$freq$conf$low)))
    out_matrix$int$up <- rbind(as.data.frame(as.matrix(object$Bayes$cred$up)),
                               as.data.frame(as.matrix(object$freq$conf$up)))
    out_matrix$omega.interval <- object$omega.interval
    out_matrix$omega.pfa <- object$freq$omega.pfa
    out_matrix$omega.error <- object$freq$omega.error
    out_matrix$omega.item.error <- object$freq$omega.item.error
    out_matrix$n.iter <- object$n.iter
    out_matrix$n.burnin <- object$n.burnin
    out_matrix$n.boot <- object$n.boot
    out_matrix$thin <- object$thin
    out_matrix$n.chains <- object$n.chains
    out_matrix$ifitem$bay_est <- object$Bayes$ifitem$est
    out_matrix$ifitem$bay_cred <- object$Bayes$ifitem$cred
    out_matrix$ifitem$freq_tab <- object$freq$ifitem
    out_matrix$para.boot <- object$para.boot
    out_matrix$inv.mat <- object$freq$inv.mat

  } else if (!is.null(object$Bayes)) {
    out_matrix$est <- as.data.frame(as.matrix(object$Bayes$est))
    out_matrix$int$low <- as.data.frame(as.matrix(object$Bayes$cred$low))
    out_matrix$int$up <- as.data.frame(as.matrix(object$Bayes$cred$up))
    out_matrix$n.iter <- object$n.iter
    out_matrix$n.burnin <- object$n.burnin
    out_matrix$thin <- object$thin
    out_matrix$n.chains <- object$n.chains
    out_matrix$ifitem$bay_est <- object$Bayes$ifitem$est
    out_matrix$ifitem$bay_cred <- object$Bayes$ifitem$cred

  } else if (!is.null(object$freq)) {
    out_matrix$est <- as.data.frame(as.matrix(object$freq$est))
    out_matrix$int$low <- as.data.frame(as.matrix(object$freq$conf$low))
    out_matrix$int$up <- as.data.frame(as.matrix(object$freq$conf$up))
    out_matrix$n.boot <- object$n.boot
    out_matrix$ifitem$freq_tab <- object$freq$ifitem
    out_matrix$omega.interval <- object$omega.interval
    out_matrix$omega.pfa <- object$freq$omega.pfa
    out_matrix$omega.error <- object$freq$omega.error
    out_matrix$omega.item.error <- object$freq$omega.item.error
    out_matrix$para.boot <- object$para.boot
    out_matrix$inv.mat <- object$freq$inv.mat

  } else {
    return(warning("no estimates calculated"))
  }
  out_matrix$call <- object$call
  out_matrix$interval <- object$interval
  out_matrix$estimates <- object$estimates
  out_matrix$complete <- object$complete
  out_matrix$pairwise <- object$miss_pairwise

  class(out_matrix) <- "summary.strel"
  out_matrix
}

#'@export
print.summary.strel <- function(x, ...) {
  n_row <- length(x$est$V1)
  mat <- data.frame(matrix(0, n_row, 0))
  mat[, 1] <- x$est
  mat[, 2] <- "   "
  mat[, 3] <- x$int$low
  mat[, 4] <- x$int$up
  row.names(mat) <- row.names(x$est)
  colnames(mat) <- c("point est", "", paste0(x$interval * 100, "% CI lower"),
                     paste0(x$interval * 100, "% CI upper"))

  cat("Call: \n")
  print.default(x$call)
  cat("\n")
  cat("Results: \n")
  print(mat, right = FALSE)
  cat("\n")
  if (!is.null(x$n.iter)) {
    cat("Bayesian point est is the posterior mean \n")
  }
  if (length(grep("freq", x$est)) > 0) {

    if (!is.null(x$inv.mat)) {
      cat("the number of bootstrap samples reduced to ")
      cat(x$inv.mat)
      cat(" because some bootstrapped matrices were singular\n")
    }
    if ("omega" %in% x$estimates) {
      if (!is.null(x$omega.pfa) & !is.null(x$omega.error)) {
        cat("frequentist omega method is pfa\n")
        cat("omega confidence interval is estimated with bootstrap because the cfa did not find a solution\n")
      }
    }

  }

  if (!is.null(x$complete)) {
    cat("Missing data handling: using listwise deletion the number of complete cases is\n")
    cat(x$complete)
  }
  if (!is.null(x$pairwise)) {
    cat("Missing data handling: using pairwise complete cases\n")
  }


  if (!is.null(x$ifitem$bay_est)) {
    n_row <- length(unlist(x$ifitem$bay_est[1])) + 1
    n_col <- 3
    names <- NULL
    for(z in 1:(n_row - 1)){
      names[z] <- paste0("x", z)
    }
    row_names <- c("original", names)

    for (i in seq_len(length(x$estimates))) {
      mat_ifitem_bay <- data.frame(matrix(0, n_row, n_col))
      row.names(mat_ifitem_bay) <- row_names

      mat_ifitem_bay[1, ] <- c(unlist(x$est)[i], unlist(x$int$low)[i], unlist(x$int$up)[i])
      mat_ifitem_bay[2:n_row, ] <- cbind(as.vector(unlist(x$ifitem$bay_est[i])),
                                         matrix(unlist(x$ifitem$bay_cred[i]), n_row - 1, 2))
      colnames(mat_ifitem_bay) <- c("point.estimate", "interval.low", "interval.up")
      cat("\n")
      cat(paste0("Bayesian ", x$estimate[i], " if item dropped: \n"))
      print(mat_ifitem_bay)

    }
  }

  if (!is.null(x$ifitem$freq_tab)) {
    n_row <- length(unlist(x$ifitem$freq_tab[1])) + 1
    n_col <- length(x$estimates)
    names <- NULL
    for (z in 1:(n_row-1)) {
      names[z] <- paste0("x", z)
    }
    row_names <- c("original", names)
    mat_ifitem_freq <- data.frame(matrix(0, n_row, n_col))
    mat_ifitem_freq[1, ] <- as.vector(unlist(x$est)[grep("freq", rownames(x$est))])
    for (i in 1:n_col) {
      mat_ifitem_freq[2:n_row, i] <- as.vector(unlist(x$ifitem$freq_tab[i]))
    }
    colnames(mat_ifitem_freq) <- x$estimates
    row.names(mat_ifitem_freq) <- c("original", names)

    cat("\n")
    cat("Frequentist point estimate if item dropped: \n")
    print(mat_ifitem_freq)

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
  out <- cbind(as.numeric(sprintf("%.5f", c(x$omega_t$mean, x$omega_t$cred))),
               as.numeric(sprintf("%.5f", c(x$omega_h$mean, x$omega_h$cred))))
  colnames(out) <- c("omega_t", "omega_h")
  row.names(out) <- c("posterior mean", paste0(x$interval * 100, "% CI lower"),
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
  out <- cbind(as.numeric(sprintf("%.5f", c(x$omega_t$est, x$omega_t$conf))),
               as.numeric(sprintf("%.5f", c(x$omega_h$est, x$omega_h$conf))))
  colnames(out) <- c("omega_t", "omega_h")
  row.names(out) <- c("point est", paste0(x$interval * 100, "% CI lower"),
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
    measures <- unname(x$model$fit.measures[c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")])
    measures <- as.numeric(sprintf("%.5f", measures))
    names(measures) <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
    cat("\nFit measures:\n")
    print.default(c(measures))
  }

}


