

test_that("Bayesian omegas are correct, missing pairwise", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2)

  expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
               c(0.8624918, 0.8459065, 0.8793787, 0.6401878, 0.5973778, 0.6870752), tolerance = 1e-3)


})

test_that("Bayesian omegas are correct, missing listwise, model sytnax specified", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  mod <- "
  f1 =~ U17_r+U22_r+U29_r+U34_r
  f2 =~ U4+U14+U19+U27
  f3 =~ U6 +U16+U28+U48
  f4 =~ U23_r +U31_r +U36_r +U46_r
  f5 =~ U10_r +U20_r +U35_r +U52_r
  "
  ee <- Bayesrel::bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2,
                          missing = "listwise", model = mod)

  expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
               c(0.8629340, 0.8432947, 0.8822594, 0.6428649, 0.5897601, 0.6868880), tolerance = 1e-3)


})


test_that("Frequentist omegas are correct with higher order model, missing pairwise", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::omegasCFA(upps, n.factors = 5)

  expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
               c(0.864759, 0.8466538, 0.8828642, 0.643837, 0.5943180, 0.6933561), tolerance = 1e-3)


})

test_that("Frequentist omegas are correct with bifactor model, missing listwise", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::omegasCFA(upps, n.factors = 5, model.type = "bi-factor", missing = "listwise")

  expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
               c(0.8718139, 0.8518070, 0.8918207, 0.6309915, 0.5792708, 0.6827122), tolerance = 1e-3)

})


test_that("Frequentist omegas are correct with bifactor model, missing listwise, model sytnax specified,
          first and last item switched", {

  data(upps, package = "Bayesrel")
  set.seed(1234)
  mod <- "
  f1 =~ U52_r+U22_r+U29_r+U34_r
  f2 =~ U4+U14+U19+U27
  f3 =~ U6 +U16+U28+U48
  f4 =~ U23_r +U31_r +U36_r +U46_r
  f5 =~ U10_r +U20_r +U35_r + U17_r
  "
  ee <- Bayesrel::omegasCFA(upps, n.factors = 5, model.type = "bi-factor", missing = "listwise",
                            model = mod)

  expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
               c(0.8702046, 0.8499908, 0.8904183, 0.6275952, 0.5760795, 0.6791109), tolerance = 1e-3)

})
