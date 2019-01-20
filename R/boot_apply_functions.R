# the basic functions for calculating and bootstrapping the internal consistency estimates

#######  measures functions ##########

applyalpha <- function(M){
  p <- ncol(M)
  a <- (p/(p-1))*(1-(sum(diag((M)))/sum(M)))
  return(a)
}

applyl2 <- function(M){
  p <- ncol(M)
  M0 <- M
  diag(M0) <- 0
  l2 <- (sum(M0) + sqrt(p/(p-1) * sum(M0^2))) / sum(M)
  return(l2)
}

applyl4 <- function(M){
  if (ncol(M) < 15) {l4 <- MaxSplitExhaustive(M)}
  else {l4 <- quant.lambda4(M)}
  return(l4)
}


applyl6 <- function(M){
  M <- cov2cor(M)
  smc <- 1 - (1 / diag(solve(M)))
  l6 <- 1 - (sum(1 - (smc)) / sum(M))
  return(l6)
}

applyglb <- function(M){
  gl <- glb.algebraic2(M)$glb
  return(gl)
}

applyomega_cfa_cov <- function(C){
  out <- omegaFreqCov(C)
  om <- out$omega
  return(om)
}

applyomega_cfa_data <- function(data){
  out <- omegaFreqData(data)
  om <- out$omega
  return(om)
}


applyomega_pa <- function(m){
  f <- princFac(m)
  l.fa <- f$loadings
  er.fa <- f$err.var
  om <- sum(l.fa)^2 / (sum(l.fa)^2 + sum(er.fa))
  if (om < 0 || om > 1 || is.na(om)) om <- NA
  return(om)
}


omegaBasic <- function(l, e){
  o <- sum(l)^2 / (sum(l)^2 + sum(e))
  return(o)
}


