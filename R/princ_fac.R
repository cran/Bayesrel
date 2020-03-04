# function to do a principal factor analysis, also principal axis method
# source:
# https://www.r-bloggers.com/iterated-principal-factor-method-of-factor-analysis-with-r/

princFac <- function(m, max_iter = 50){
  # r <- cov2cor(m)
  r <- m
  r_smc <- (1 - 1 / diag(try(solve(cov2cor(r)))))
  diag(r) <- r_smc
  min_error <- .001
  com_iter <- c()
  h2 <- sum(diag(r))
  error <- h2
  i <- 1
  while (error > min_error || i == 1){
    r_eigen <- eigen(r)

    lambda <- as.matrix(r_eigen$vectors[, 1] * sqrt(r_eigen$values[1]))

    r_mod <- lambda %*% t(lambda)
    r_mod_diag <- diag(r_mod)

    h2_new <- sum(r_mod_diag)
    error <- abs(h2 - h2_new)

    h2 <- h2_new

    com_iter <- append(com_iter, h2_new)
    diag(r) <- r_mod_diag
    i <- i + 1
    if (i > max_iter) {
      error <- 0
    }
  }

  # h2 <- rowSums(lambda^2)
  # u2 <- 1 - h2
  # com <- rowSums(lambda^2)^2 / rowSums(lambda^4)
  # pf.loadings <- data.frame(cbind(round(lambda, 3), round(h2, 3), round(u2, 3), round(com, 3)))
  # colnames(pf.loadings) <- c("factor", "h2", "u2", "com")
  if(sum(lambda) < 0){
    lambda <- -lambda
  }
  L <- lambda %*% t(lambda)
  e <- diag(m - L)
  return(list(loadings = lambda, err_var = e))
}

# source:
# from the psych package:
# Revelle, W. (2018) psych: Procedures for Personality and Psychological Research, Northwestern University, Evanston,
# Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.8.4.
corSmooth2 <- function (x, eig_tol = 10^-12) {
  eigens <- try(eigen(x), TRUE)
  if (inherits(eigens, "try-error")) {
    warning("I am sorry, there is something seriously wrong with the correlation matrix,\ncor_smooth failed to smooth it because some of the eigen values are NA.  \nAre you sure you specified the data correctly?")
  }
  else {
    if (min(eigens$values) < .Machine$double.eps) {
      warning("Matrix was not positive definite, smoothing was done")
      eigens$values[eigens$values < eig_tol] <- 100 *
        eig_tol
      nvar <- dim(x)[1]
      tot <- sum(eigens$values)
      eigens$values <- eigens$values * nvar/tot
      cnames <- colnames(x)
      rnames <- rownames(x)
      x <- eigens$vectors %*% diag(eigens$values) %*%
        t(eigens$vectors)
      x <- cov2cor(x)
      colnames(x) <- cnames
      rownames(x) <- rnames
    }
  }
  return(x)
}
