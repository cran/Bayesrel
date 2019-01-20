

# this function calls on other functions in order to return the frequentist estimates
# and parametric bootstrapped confidence intervals, sampling from a multivariate normal distribution

freqFun_para <- function(data, boot.n, estimates, interval, omega.freq.method,
                         item.dropped, alpha.int.analytic, omega.fit){
  p <- ncol(data)
  n <- nrow(data)
  res <- list()
  res$covsamp <- NULL
  if ("alpha" %in% estimates || "lambda2" %in% estimates || "lambda4" %in% estimates || "lambda6" %in% estimates ||
      "glb" %in% estimates || omega.freq.method == "pa"){
    boot.data <- array(0, c(boot.n, n, p))
    boot.cov <- array(0, c(boot.n, p, p))
    for (i in 1:boot.n){
      boot.data[i, , ] <- MASS::mvrnorm(n, colMeans(data), cov(data))
      boot.cov[i, , ] <- cov(boot.data[i, , ])
    }
    res$covsamp <- boot.cov
  }
  if (item.dropped){
    Ctmp <- array(0, c(p, p - 1, p - 1))
    Dtmp <- array(0, c(p, n, p - 1))
    for (i in 1:p){
      Ctmp[i, , ] <- cov(data)[-i, -i]
      Dtmp[i, , ] <- data[, -i]
    }
  }
  if ("alpha" %in% estimates){
    res$est$freq.alpha <- applyalpha(cov(data))
    if (alpha.int.analytic){
      int <- ciAlpha(1 - interval, n, cov(data))
      res$conf$low$freq.alpha <- int[1]
      res$conf$up$freq.alpha <- int[2]

    } else{
      alpha.obj <- apply(boot.cov, 1, applyalpha)
      if (length(unique(round(alpha.obj, 4))) == 1){
        res$conf$low$freq.alpha <- 1
        res$conf$up$freq.alpha <- 1
      }
      else{
        res$conf$low$freq.alpha <- quantile(alpha.obj, probs = (1 - interval)/2, na.rm = T)
        res$conf$up$freq.alpha <- quantile(alpha.obj, probs = interval + (1 - interval)/2, na.rm = T)
      }
      res$boot$alpha <- alpha.obj
    }
    if (item.dropped){
      res$ifitem$alpha <- apply(Ctmp, 1, applyalpha)
    }
  }
  if ("lambda2" %in% estimates){
    res$est$freq.l2 <- applyl2(cov(data))
    l2.obj <- apply(boot.cov, 1, applyl2)
    if (length(unique(round(l2.obj, 4))) == 1){
      res$conf$low$freq.l2 <- 1
      res$conf$up$freq.l2 <- 1
    }
    else{
      res$conf$low$freq.l2 <- quantile(l2.obj, probs = (1 - interval)/2, na.rm = T)
      res$conf$up$freq.l2 <- quantile(l2.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$l2 <- l2.obj
    if (item.dropped){
      res$ifitem$l2 <- apply(Ctmp, 1, applyl2)
    }
  }

  if ("lambda4" %in% estimates){
    res$est$freq.l4 <- applyl4(cov(data))
    l4.obj <- apply(boot.cov, 1, applyl4)
    if (length(unique(round(l4.obj, 4))) == 1){
      res$conf$low$freq.l4 <- 1
      res$conf$up$freq.l4 <- 1
    }
    else{
      res$conf$low$freq.l4 <- quantile(l4.obj, probs = (1 - interval)/2, na.rm = T)
      res$conf$up$freq.l4 <- quantile(l4.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$l4 <- l4.obj
    if (item.dropped){
      res$ifitem$l4 <- apply(Ctmp, 1, applyl4)
    }
  }

  if ("lambda6" %in% estimates){
    res$est$freq.l6 <- applyl6(cov(data))
    l6.obj <- apply(boot.cov, 1, applyl6)
    if (length(unique(round(l6.obj, 4))) == 1){
      res$conf$low$freq.l6 <- 1
      res$conf$up$freq.l6 <- 1
    }
    else{
      res$conf$low$freq.l6 <- quantile(l6.obj, probs = (1 - interval)/2, na.rm = T)
      res$conf$up$freq.l6 <- quantile(l6.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$l6 <- l6.obj
    if (item.dropped){
      res$ifitem$l6 <- apply(Ctmp, 1, applyl6)
    }
  }
  if ("glb" %in% estimates){
    res$est$freq.glb <- applyglb(cov(data))
    glb.obj <- apply(boot.cov, 1, applyglb)
    if (length(unique(round(glb.obj, 4))) == 1){
      res$conf$low$freq.glb <- 1
      res$conf$up$freq.glb <- 1
    }
    else{
      res$conf$low$freq.glb <- quantile(glb.obj, probs = (1 - interval)/2, na.rm = T)
      res$conf$up$freq.glb <- quantile(glb.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$glb <- glb.obj
    if (item.dropped){
      res$ifitem$glb <- apply(Ctmp, 1, applyglb)
    }
  }

  #omega --------------------------------------------------------------------------
  if ("omega" %in% estimates){
    if (omega.freq.method == "cfa"){
      out <- omegaFreqData(data)
      res$est$freq.omega <- out$omega
      res$loadings <- out$loadings
      res$resid.var <- out$errors
      res$conf$low$freq.omega <- out$omega.low
      res$conf$up$freq.omega <- out$omega.up

      if (omega.fit) {res$fit$omega <- out$indices}
      if (item.dropped){
        res$ifitem$omega <- apply(Dtmp, 1, applyomega_cfa_data)
      }
    }
    if (omega.freq.method == "pa"){
      res$est$freq.omega <- applyomega_pa(cov(data))
      omega.obj <- apply(boot.cov, 1, applyomega_pa)
      if (length(unique(round(omega.obj, 4))) == 1){
        res$conf$low$freq.omega <- 1
        res$conf$up$freq.omega <- 1
      }
      else{
        res$conf$low$freq.omega <- quantile(omega.obj, probs = (1 - interval)/2, na.rm = T)
        res$conf$up$freq.omega <- quantile(omega.obj, probs = interval + (1 - interval)/2, na.rm = T)
      }
      res$boot$omega <- omega.obj
      if (item.dropped){
        res$ifitem$omega <- apply(Ctmp, 1, applyomega_pa)
      }
    }
  }
  return(res)
}
