# function to do a principal factor analysis, also principal axis method
# source:
# https://www.r-bloggers.com/iterated-principal-factor-method-of-factor-analysis-with-r/

princFac <- function(m, max.iter = 50){
  r <- cov2cor(m)
  r.smc <- (1 - 1 / diag(solve(r)))
  diag(r) <- r.smc
  min.error <- .001
  com.iter <- c()
  h2 <- sum(diag(r))
  error <- h2
  i <- 1
  while (error > min.error || i == 1){
    r.eigen <- eigen(r)

    lambda <- as.matrix(r.eigen$vectors[, 1] * sqrt(r.eigen$values[1]))

    r.mod <- lambda %*% t(lambda)
    r.mod.diag <- diag(r.mod)

    h2.new <- sum(r.mod.diag)
    error <- abs(h2 -h2.new)

    h2 <- h2.new

    com.iter <- append(com.iter, h2.new)
    diag(r) <- r.mod.diag
    i <- i + 1
    if (i > max.iter) {
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
  e <- 1 - lambda^2
  return(list(loadings = lambda, err.var = e))
}

# source:
# from the psych package:
# Revelle, W. (2018) psych: Procedures for Personality and Psychological Research, Northwestern University, Evanston,
# Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.8.4.
corSmooth2 <- function (x, eig.tol = 10^-12) {
  eigens <- try(eigen(x), TRUE)
  if (inherits(eigens, "try-error")) {
    warning("I am sorry, there is something seriously wrong with the correlation matrix,\ncor.smooth failed to smooth it because some of the eigen values are NA.  \nAre you sure you specified the data correctly?")
  }
  else {
    if (min(eigens$values) < .Machine$double.eps) {
      warning("Matrix was not positive definite, smoothing was done")
      eigens$values[eigens$values < eig.tol] <- 100 *
        eig.tol
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
