


set.seed(1234)

test_that("Estimates lambda2 and omega are correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::strel(asrm, estimates = c("lambda2", "omega"), n.iter = 500, n.boot = 200)

  expect_equal(ee$Bayes$est$Bayes_lambda2, 0.7970929, tolerance = 1e-3)
  expect_equal(ee$freq$est$freq_lambda2, 0.7960336, tolerance = 1e-3)
  expect_equal(ee$Bayes$est$Bayes_omega, 0.7735565, tolerance = 1e-3)
  expect_equal(ee$freq$est$freq_omega, 0.7919616, tolerance = 1e-3)

  expect_equal(ee$Bayes$cred$low$Bayes_omega, 0.6929867, tolerance = 1e-3)
  if (as.numeric(R.Version()$minor) > 6) {
    expect_equal(as.numeric(ee$freq$conf$up$freq_lambda2), 0.867655, tolerance = 1e-3)
  } # because of the change in the RNG brought by the new R version

})



test_that("Bayes glb is correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::strel(asrm, estimates = "glb", n.iter = 500, freq = F)

  expect_equal(ee$Bayes$est$Bayes_glb, 0.8589229, tolerance = 1e-3)
  expect_equal(ee$Bayes$cred$up$Bayes_glb, 0.9037633, tolerance = 1e-3)


})


test_that("Bayes Alpha if item deleted is correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::strel(asrm, estimates = "alpha", n.iter = 500, freq = F, item.dropped = T)
  expect_equal(as.numeric(ee$Bayes$ifitem$est$alpha), c(0.7225611, 0.7311794, 0.7935113, 0.7658315, 0.7323329),
               tolerance = 1e-3)

})




test_that("Bayes prior and posterior probability for Alpha >.8 is correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  tt <- Bayesrel::strel(asrm, estimates = "alpha", n.iter = 500, freq = F)
  ee <- Bayesrel::p_strel(tt, "alpha", .8)

  expect_equal(as.numeric(ee), c(0.1552618, 0.3711111), tolerance = 1e-3)

})
