# this function calls on other functions in order to return the sampled estimates
# and the credible intervals together with the posterior distribution objects
# to be passed on for forther analysis

gibbsFun <- function(data, n.iter, n.burnin, estimates, interval, item.dropped, omega.fit){
  p <- ncol(data)
  res <- list()
  if ("alpha" %in% estimates || "lambda2" %in% estimates || "lambda4" %in% estimates || "lambda6" %in% estimates ||
      "glb" %in% estimates){
    C <- covSamp(data, n.iter, n.burnin)
    # Cmed <- apply(C, c(2, 3), median)
    if (item.dropped) {
      Ctmp <- array(0, c(p, n.iter - n.burnin, p - 1, p - 1))
      for (i in 1:p){
        Ctmp[i, , , ] <- C[, -i, -i]
      }
    }
  } else {
    C = NULL
  }
  res$covsamp <- C

  if ("alpha" %in% estimates){
    res$samp$bayes.alpha <- coda::as.mcmc(apply(C, MARGIN = 1, applyalpha))
    int <- coda::HPDinterval(res$samp$bayes.alpha, prob = interval)
    res$cred$low$bayes.alpha <- int[1]
    res$cred$up$bayes.alpha <- int[2]
    res$est$bayes.alpha<- median(res$samp$bayes.alpha)
    # res$est$bayes.alpha<- applyalpha(Cmed)
    if (item.dropped){
      res$ifitem$samp$alpha <- apply(Ctmp, c(1, 2), applyalpha)
      res$ifitem$est$alpha <- apply(res$ifitem$samp$alpha, 1, median)
    }
  }

  if ("lambda2" %in% estimates){
    res$samp$bayes.l2 <- coda::as.mcmc(apply(C, MARGIN = 1, applyl2))
    int <- coda::HPDinterval(res$samp$bayes.l2, prob = interval)
    res$cred$low$bayes.l2 <- int[1]
    res$cred$up$bayes.l2 <- int[2]
    res$est$bayes.l2<- median(res$samp$bayes.l2)
    # res$est$bayes.l2<- applyl2(Cmed)
    if (item.dropped){
      res$ifitem$samp$l2 <- apply(Ctmp, c(1, 2), applyl2)
      res$ifitem$est$l2 <- apply(res$ifitem$samp$l2, 1, median)
    }
  }

  if ("lambda4" %in% estimates){
    res$samp$bayes.l4 <- coda::as.mcmc(apply(C, MARGIN = 1, applyl4))
    int <- coda::HPDinterval(res$samp$bayes.l4, prob = interval)
    res$cred$low$bayes.l4 <- int[1]
    res$cred$up$bayes.l4 <- int[2]
    res$est$bayes.l4<- median(res$samp$bayes.l4)
    # res$est$bayes.l4<- applyl4(Cmed)
    if (item.dropped){
      res$ifitem$samp$l4 <- apply(Ctmp, c(1, 2), applyl4)
      res$ifitem$est$l4 <- apply(res$ifitem$samp$l4, 1, median)
    }
  }

  if ("lambda6" %in% estimates){
    res$samp$bayes.l6 <- coda::as.mcmc(apply(C, MARGIN = 1, applyl6))
    int <- coda::HPDinterval(res$samp$bayes.l6, prob = interval)
    res$cred$low$bayes.l6 <- int[1]
    res$cred$up$bayes.l6 <- int[2]
    res$est$bayes.l6<- median(res$samp$bayes.l6)
    # res$est$bayes.l6<- applyl6(Cmed)
    if (item.dropped){
      res$ifitem$samp$l6 <- apply(Ctmp, c(1, 2), applyl6)
      res$ifitem$est$l6 <- apply(res$ifitem$samp$l6, 1, median)
    }
  }

  if ("glb" %in% estimates){
    res$samp$bayes.glb <- coda::as.mcmc(apply(C, MARGIN = 1, applyglb))
    int <- coda::HPDinterval(res$samp$bayes.glb, prob = interval)
    res$cred$low$bayes.glb <- int[1]
    res$cred$up$bayes.glb <- int[2]
    res$est$bayes.glb<- median(res$samp$bayes.glb)
    # res$est$bayes.glb<- applyglb(Cmed)
    if (item.dropped){
      res$ifitem$samp$glb <- apply(Ctmp, c(1, 2), applyglb)
      res$ifitem$est$glb <- apply(res$ifitem$samp$glb, 1, median)
    }
  }

  # special case omega -----------------------------------------------------------------
  if ("omega" %in% estimates){
    om.samp <- omegaSampler(data, n.iter, n.burnin)
    res$samp$bayes.omega <- coda::as.mcmc(om.samp$omega)
    res$loadings <- apply(om.samp$lambda, 2, median)
    res$resid.var <- apply(om.samp$psi, 2, median)
    int <- coda::HPDinterval(res$samp$bayes.omega, prob = interval)
    res$cred$low$bayes.omega <- int[1]
    res$cred$up$bayes.omega<- int[2]
    res$est$bayes.omega <- median(res$samp$bayes.omega)
    if (omega.fit){
      ppcOmega(data, res$loadings, res$resid.va)
      res$fit$omega <- recordPlot()
      # plot.new()
    }
    if (item.dropped){
      om.samp.ifitem <- matrix(0, p, n.iter - n.burnin)
      for (i in 1:p){
        tmp <- data[-i, -i]
        om.samp.ifitem[i, ] <- omegaSampler(tmp, n.iter, n.burnin)$omega
      }
      res$ifitem$samp$omega <- om.samp.ifitem
      res$ifitem$est$omega <- apply(om.samp.ifitem, 1, mean)
    }
  }

  return(res)

}
