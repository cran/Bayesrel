
tol <- ifelse(Sys.getenv("HOME") == "/Users/julius" || Sys.getenv("IAMHOME") == TRUE, 1e-5, 1e-2)

data(asrm, package = "Bayesrel")

# Estimates lambda2 and omega are correct
set.seed(1234)
ee <- strel(asrm, estimates = c("lambda2", "omega"), n.iter = 100, n.boot = 100, n.chains = 2,
                      disableMcmcCheck = TRUE)

expect_equal(ee$Bayes$est$Bayes_lambda2, 0.7948049, tolerance = tol)
expect_equal(ee$freq$est$freq_lambda2, 0.7960336, tolerance = tol)
expect_equal(ee$Bayes$est$Bayes_omega, 0.7708523, tolerance = tol)
expect_equal(ee$freq$est$freq_omega, 0.7919616, tolerance = tol)
expect_equal(ee$Bayes$cred$low$Bayes_omega, 0.6719616, tolerance = tol)
expect_equal(ee$freq$conf$low$freq_omega, 0.7194611, tolerance = tol)
expect_equal(as.numeric(ee$freq$conf$up$freq_lambda2), 0.865121, tolerance = tol)



# Bayes glb is correct
set.seed(1234)
ee <- strel(asrm, estimates = "glb", n.iter = 100, freq = FALSE, n.chains = 1,
                      disableMcmcCheck = TRUE)

expect_equal(ee$Bayes$est$Bayes_glb, 0.8542316, tolerance = tol)
expect_equal(ee$Bayes$cred$up$Bayes_glb, 0.8950283, tolerance = tol)



# Bayes alpha and omega if item deleted is correct
set.seed(1234)
ee <- strel(asrm, estimates = c("alpha", "omega"), n.iter = 100, freq = FALSE, item.dropped = TRUE,
                      n.chains = 2, disableMcmcCheck = TRUE)
expect_equal(as.numeric(ee$Bayes$ifitem$est$alpha[1:2]), c(0.7207363, 0.7245768),
             tolerance = tol)
expect_equal(as.numeric(ee$Bayes$ifitem$cred$alpha[c(1, 10)]), c(0.6450673, 0.8049170),
             tolerance = tol)
expect_equal(as.numeric(ee$Bayes$ifitem$est$omega[4:5]), c(0.7449255, 0.7204507),
             tolerance = tol)
expect_equal(as.numeric(ee$Bayes$ifitem$cred$omega[c(1, 10)]), c(0.6009078, 0.8274682),
             tolerance = tol)



# Freq omega with PFA is correct
set.seed(1234)
tt <- strel(asrm, estimates = "omega", n.boot = 100, Bayes = FALSE, omega.freq.method = "pfa")
expect_equal(as.numeric(tt$freq$est$freq_omega), c(0.7966209), tolerance = tol)
if (as.numeric(R.Version()$major >= 4)) {
  expect_equal(as.numeric(tt$freq$conf$up$freq_omega), 0.8617495, tolerance = tol)
} # because of the change in the RNG brought by the new R version



# Bayes prior and posterior probability for Alpha >.8 is correct
set.seed(1234)
ee <- strel(asrm, estimates = "alpha", n.iter = 100, freq = FALSE, n.chains = 2,
                      disableMcmcCheck = TRUE)
tt <- pStrel(ee, "alpha", .8)

expect_equal(as.numeric(tt), c(0.1556125, 0.3300000), tolerance = tol)


# Omega results with missing data are correct
set.seed(1234)
ee <- strel(asrm_mis, estimates = c("omega"), n.iter = 100, n.chains = 2, n.boot = 100,
                      disableMcmcCheck = TRUE)
expect_equal(as.numeric(ee$Bayes$cred$low$Bayes_omega), c(0.6900274),
             tolerance = tol)
expect_equal(as.numeric(ee$freq$est$freq_omega), c(0.7943602),
             tolerance = tol)


# Frequentist Lambda6 results with missing data and parametric bootstrap are correct
set.seed(1234)
ee <- strel(asrm_mis, estimates = c("lambda6"), Bayes = FALSE, n.boot = 100, para.boot = TRUE)
expect_equal(as.numeric(ee$freq$est$freq_lambda6), c(0.7927271),
             tolerance = tol)
expect_equal(as.numeric(ee$freq$conf$low$freq_lambda6), 0.7188984, tolerance = tol)


# Results with input cov matrix are correct
cc <- cov(asrm)
set.seed(1234)
ee <- strel(cov.mat = cc, estimates = c("lambda2"), n.iter = 100, n.chains = 2, n.boot = 100, n.obs = 500,
                      disableMcmcCheck = TRUE)
expect_equal(as.numeric(ee$Bayes$cred$up$Bayes_lambda2), c(0.8215358),
             tolerance = tol)
expect_equal(as.numeric(ee$freq$est$freq_lambda2), c(0.7960336),
             tolerance = tol)
if (Sys.info()["sysname"] == "Linux") tol <- 1e-2
expect_equal(as.numeric(ee$freq$conf$low$freq_lambda2), 0.771300720790473, tolerance = tol)
tol <- 1e-3

# Frequentist glb is correct
set.seed(1234)
ee <- strel(asrm, estimates = "glb", n.boot = 100, freq = TRUE, Bayes = FALSE)

expect_equal(c(ee$freq$est$freq_glb), 0.8459531, tolerance = tol)
expect_equal(c(ee$freq$conf$up$freq_glb, use.names = FALSE), 0.9011926, tolerance = tol)



# Bayesian estimates lambda2 and omega are correct with adjusted priors
set.seed(1234)
ee <- strel(asrm_mis, estimates = c("lambda2", "omega"), n.iter = 100, freq = FALSE, n.chains = 2,
                      k0 = 1, df0 = 10, a0 = 6, b0 = 10, m0 = 1, missing = "listwise", item.dropped = TRUE,
                      disableMcmcCheck = TRUE)

expect_equal(ee$Bayes$est$Bayes_lambda2, 0.7892705, tolerance = tol)
expect_equal(ee$Bayes$est$Bayes_omega, 0.6815307, tolerance = tol)
expect_equal(ee$Bayes$cred$low$Bayes_omega, 0.543211, tolerance = tol)
expect_equal(ee$Bayes$cred$up$Bayes_lambda2, 0.8596607, tolerance = tol)

expect_equal(ee$Bayes$ifitem$est$lambda2, c(0.7134186, 0.7379066, 0.7846232, 0.7670343, 0.7335534), tolerance = tol)
expect_equal(ee$Bayes$ifitem$est$omega, c(0.5280521, 0.6051296, 0.6450282, 0.5913725, 0.5508041), tolerance = tol)
expect_equal(c(ee$Bayes$ifitem$cred$lambda2), c(0.6296647, 0.6475465, 0.7102916, 0.6876639, 0.6249340, 0.8008962,
                                                0.8211305, 0.8507353, 0.8494441, 0.8130506), tolerance = tol)
expect_equal(c(ee$Bayes$ifitem$cred$omega), c(0.3427576, 0.4209486, 0.5482522, 0.3747722, 0.3586816, 0.7395393,
                                              0.7367284, 0.7716797, 0.7647747, 0.7282343), tolerance = tol)

set.seed(1234)
tt <- pStrel(ee, "lambda2", .7)
expect_equal(as.numeric(tt), c(0.1052664, 0.99000000), tolerance = tol)

tt2 <- pStrel(ee, "omega", .7)
expect_equal(as.numeric(tt2), c(0.5186924, 0.4600000), tolerance = tol)



# Omega fit is correct
set.seed(1234)
ee <- strel(asrm, estimates = "omega", n.iter = 200, n.chains = 2)
ff <- omegaFit(ee, asrm, ppc = FALSE, cutoff = .05)
expect_equal(unlist(ff, use.names = FALSE), c(13.06419359, 0.06170546, 0.13031834, 0.06918759, 0.22394139, 0.05666667,
                                            12.02439250, 5.00000000, 0.03445507, 0.13420605, 0.03365073, 0.23333645,
                                            0.05875407),
             tolerance = tol)
