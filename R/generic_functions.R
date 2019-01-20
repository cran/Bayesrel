
#'@export
print.strel <- function(x, ...){
  if (!is.null(x$freq)){
    est <- cbind(as.data.frame(as.matrix(x$bay$est)),
                 as.data.frame(as.matrix(x$freq$est)))
    colnames(est) <- c("bayes", "frequentist")
  }
  else{
    est <- as.data.frame(as.matrix(x$bay$est))
    colnames(est) <- "bayes"
  }
  row.names(est) <- x$estimates
  cat("Call: \n")
  print.default(x$call)
  cat("\n")
  cat("Point Estimates of Internal Consistency Measures: \n")
  cat("\n")
  print(est)

}

#'@export
summary.strel <- function(object, ...){

  out.matrix <- list()
  if (!is.null(object$freq)){
    out.matrix$est <- rbind(as.data.frame(as.matrix(object$bay$est)),
                            as.data.frame(as.matrix(object$freq$est)))
    out.matrix$int$low <- rbind(as.data.frame(as.matrix(object$bay$cred$low)),
                                as.data.frame(as.matrix(object$freq$conf$low)))
    out.matrix$int$up <- rbind(as.data.frame(as.matrix(object$bay$cred$up)),
                               as.data.frame(as.matrix(object$freq$conf$up)))
    out.matrix$omega.freq.method <- object$omega.freq.method
  } else{
    out.matrix$est <- as.data.frame(as.matrix(object$bay$est))
    out.matrix$int$low <- as.data.frame(as.matrix(object$bay$cred$low))
    out.matrix$int$up <- as.data.frame(as.matrix(object$bay$cred$up))
  }
  out.matrix$call <- object$call
  out.matrix$n.iter <- object$n.iter
  out.matrix$n.burnin <- object$n.burnin
  out.matrix$interval <- object$interval
  out.matrix$estimates <- object$estimates
  out.matrix$fit.indices <- object$freq$fit$omega
  out.matrix$ifitem$bay.tab <- object$bay$ifitem$est
  out.matrix$ifitem$freq.tab <- object$freq$ifitem

  class(out.matrix) <- "summary.strel"
  out.matrix
}

#'@export
print.summary.strel <- function(x, ...){
  n.row <- length(x$est$V1)
  mat <- data.frame(matrix(0, n.row, 0))
  mat[, 1] <- x$est
  mat[, 2] <- '   '
  mat[, 3] <- x$int$low
  mat[, 4] <- x$int$up
  row.names(mat) <- row.names(x$est)
  colnames(mat) <- c("estimate", '', "interval.low", "interval.up")

  cat("Call: \n")
  print.default(x$call)
  cat("Results: \n")
  print(mat, right = F)
  cat("\n")
  cat("iterations: ")
  print.default(x$n.iter)
  cat("burnin: ")
  print.default(x$n.burnin)
  cat("uncertainty interval:")
  print.default(x$interval)

  if (length(grep("freq", x$est)) > 0){
    if ("omega" %in% x$estimates){
      cat("frequentist omega method is:")
      print.default(x$omega.freq.method)
      cat("omega confidence interval is estimated with:")
      if (x$omega.freq.method == "pa") {print.default("bootstrap")}
      if (x$omega.freq.method == "cfa") {print.default("maximum likelihood z-value")}
    }
    if (!is.null(x$fit.indices)){
      options(scipen = 999)
      cat("\nFrequentist fit of 1-factor-model for omega is:\n")
      print.default(as.matrix(x$fit.indices))
    }
  }


  if (!is.null(x$ifitem$bay.tab)){
    n.row <- length(unlist(x$ifitem$bay.tab[1])) + 1
    n.col <- length(x$ifitem$bay.tab)
    mat.ifitem.bay <- data.frame(matrix(0, n.row, n.col))
    mat.ifitem.bay[1, ] <- as.vector(unlist(x$est))[1:n.col]
    for (i in 1:n.col){
      mat.ifitem.bay[2:n.row, i] <- as.vector(unlist(x$ifitem$bay.tab[i]))
    }
    colnames(mat.ifitem.bay) <- x$estimates
    names <- NULL
    for(i in 1:(n.row-1)){
      names[i] <- paste0("x", i)
    }
    row.names(mat.ifitem.bay) <- c("original", names)

    if (length(grep("freq", x$est)) > 0){
      mat.ifitem.freq <- data.frame(matrix(0, n.row, n.col))
      mat.ifitem.freq[1, ] <- as.vector(unlist(x$est)[(n.col+1):(n.col*2)])
      for (i in 1:n.col){
        mat.ifitem.freq[2:n.row, i] <- as.vector(unlist(x$ifitem$freq.tab[i]))
      }
      colnames(mat.ifitem.freq) <- x$estimates
      row.names(mat.ifitem.freq) <- c("original", names)
    }

    cat("\n")
    cat("Bayesian coefficient if item dropped: \n")
    print(mat.ifitem.bay)
    if (length(grep("freq", x$est)) > 0){
      cat("\n")
      cat("Frequentist coefficient if item dropped: \n")
      print(mat.ifitem.freq)
    }
  }

}

