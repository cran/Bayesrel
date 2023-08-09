
tol <- ifelse(Sys.getenv("HOME") == "/Users/julius" || Sys.getenv("IAMHOME") == TRUE, 1e-5, 1e-2)
data(upps, package = "Bayesrel")

# Bayesian omegas are correct, missing pairwise, and fit indices are good
set.seed(1234)
ee <- bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2, model.type = "second-order")

expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
             c(0.8633604, 0.8448201, 0.8790073, 0.6399786, 0.5956023, 0.6922816),
             tolerance = tol)

ff <- multiFit(ee, upps, ppc = FALSE, cutoff = .06)
expect_equal(as.numeric(c(ff$LR, ff$srmr_pointEst, ff$rmsea_pointEst, ff$rmsea_ci, ff$p_rmsea)),
             c(409.03892661, 0.07028219, 0.05616901, 0.05424941, 0.05760252,
               1.00000000), tolerance = tol)



# Bayesian omegas are correct, missing listwise, model sytnax specified
set.seed(1234)
mod <- "
f1 =~ U17_r+U22_r+U29_r+U34_r
f2 =~ U4+U14+U19+U27
f3 =~ U6 +U16+U28+U48
f4 =~ U23_r +U31_r +U36_r +U46_r
f5 =~ U10_r +U20_r +U35_r +U52_r
"
ee <- bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2,
                        missing = "listwise", model = mod, model.type = "second-order")

expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
             c(0.8624262, 0.8451732, 0.8801895, 0.6402192, 0.5909983, 0.6860989), tolerance = tol)


# Frequentist omegas are correct with higher order model, missing pairwise
set.seed(1234)
ee <- omegasCFA(upps, n.factors = 5)

expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
             c(0.864759, 0.8466538, 0.8828642, 0.643837, 0.5943180, 0.6933561), tolerance = tol)


# Frequentist omegas are correct with bifactor model, missing listwise

data(upps, package = "Bayesrel")
set.seed(1234)
ee <- omegasCFA(upps, n.factors = 5, model.type = "bi-factor", missing = "listwise")

expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
             c(0.8718139, 0.8518070, 0.8918207, 0.6309915, 0.5792708, 0.6827122), tolerance = tol)


# Frequentist omegas are correct with bifactor model, missing listwise, model sytnax specified,
# first and last item switched
set.seed(1234)
mod <- "
f1 =~ U52_r+U22_r+U29_r+U34_r
f2 =~ U4+U14+U19+U27
f3 =~ U6 +U16+U28+U48
f4 =~ U23_r +U31_r +U36_r +U46_r
f5 =~ U10_r +U20_r +U35_r + U17_r
"
ee <- omegasCFA(upps, n.factors = 5, model.type = "bi-factor", missing = "listwise",
                          model = mod)

expect_equal(c(ee$omega_t$est, ee$omega_t$conf, ee$omega_h$est, ee$omega_h$conf),
             c(0.8702046, 0.8499908, 0.8904183, 0.6275952, 0.5760795, 0.6791109), tolerance = tol)



# Bayesian omegas are correct with altered prior hyperparameters
set.seed(1234)
ee <- bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2,
                        a0 = 6, b0 = 10, l0 = 1, c0 = 10, d0 = 6, beta0 = 2, p0 = 12, R0 = 5,
                        model.type = "second-order")

expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
             c(0.8476525, 0.8275225, 0.8644323, 0.6255007, 0.5799303, 0.6704586),
             tolerance = tol)


# Bayesian omegas are correct with cross loadings, checks omega fit, and prior/post prob
mod <- "
f1 =~ U17_r + U22_r + U29_r + U34_r + U19
f2 =~ U4 + U14 + U19 + U27
f3 =~ U6 + U16 + U28 + U48 + U29_r
f4 =~ U23_r + U31_r + U36_r + U46_r + U34_r
f5 =~ U10_r + U20_r + U35_r + U52_r + U19
"
set.seed(1234)
ee <- bomegas(upps, n.iter = 200, n.burnin = 50, n.chains = 2,
                        model = mod, missing = "listwise", param.out = TRUE, model.type = "second-order")
ll <- apply(ee$model$lambda, c(3, 4), mean)

expect_equal(ll,
             matrix(c(-0.5236519, -0.5287243, -0.4082028, -0.628011, 0, 0, 0.07909764, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.003327511, 0.004258746, -0.0007414194,
                      0.001434998, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.001336209, 0, 0, 0, 0,
                      0, -0.0008860797, 0.003827371, 0.002945267, -0.002850434, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0.005032167, 0, 0, 0, 0, 0, 0, 0, 0, -0.004571036, 0.001070471, -0.002486544,
                      0.003125839, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.002564826, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0.003184005, 0.003274632, 0.003072611, -0.004328613), 20, 5),
             tolerance = tol)

ff <- multiFit(ee, upps, ppc = FALSE, cutoff = .06)
expect_equal(as.numeric(c(ff$LR, ff$srmr_pointEst, ff$rmsea_pointEst, ff$rmsea_ci, ff$p_rmsea)),
             c(354.90075890, 0.06332863, 0.05226046, 0.05029409, 0.05479048,
               1.00000000), tolerance = tol)

tt <- pOmegas(ee)
expect_equal(as.numeric(tt),
             c(0.2350000, 1.0000000, 0.2020000, 0.9333333), tolerance = tol)


# Frequentist omegas are correct with cross loadings

mod <- "
f1 =~ U17_r + U22_r + U29_r + U34_r + U19
f2 =~ U4 + U14 + U19 + U27
f3 =~ U6 + U16 + U28 + U48 + U29_r
f4 =~ U23_r + U31_r + U36_r + U46_r + U34_r
f5 =~ U10_r + U20_r + U35_r + U52_r + U19
"
ee <- omegasCFA(upps, model = mod, missing = "listwise")
expect_equal(as.numeric(c(unlist(ee[[1]]), unlist(ee[[2]]))),
             c(0.8648853, 0.8466552, 0.8831155, 0.6464642, 0.5962111, 0.6967174),
             tolerance = tol)


# Bayesian omegas for bi-factor model are correct, missing pairwise, fit indices are good,
# prior and posterior prob work
set.seed(1234)
ee <- bomegas(upps, n.factors = 5, n.iter = 200, n.burnin = 50, n.chains = 2, model.type = "bi-factor")

expect_equal(c(ee$omega_t$mean, ee$omega_t$cred, ee$omega_h$mean, ee$omega_h$cred),
             c(0.8677726, 0.8503604, 0.8845738, 0.6249922, 0.5724003, 0.6760027),
             tolerance = tol)

ff <- multiFit(ee, upps, ppc = FALSE, cutoff = .06)
expect_equal(as.numeric(c(ff$LR, ff$srmr_pointEst, ff$rmsea_pointEst, ff$rmsea_ci, ff$p_rmsea)),
             c(316.24190905, 0.05806143, 0.04942153, 0.04698159, 0.05242331, 1.00000000), tolerance = tol)

tt <- pOmegas(ee)
expect_equal(as.numeric(tt),
             c(0.2730, 1.0000, 0.2725, 0.8200), tolerance = tol)


# Bayesian omegas are correct from correlated model, missing data, and crossloadings
mod <- "
f1 =~ U17_r + U22_r + U29_r + U34_r + U19
f2 =~ U4 + U14 + U19 + U27
f3 =~ U6 + U16 + U28 + U48 + U29_r
f4 =~ U23_r + U31_r + U36_r + U46_r + U34_r
f5 =~ U10_r + U20_r + U35_r + U52_r + U19
"
set.seed(1234)
ee <- bomegas(upps, n.iter = 200, n.burnin = 50, n.chains = 2,
              model = mod, missing = "impute", param.out = TRUE, model.type = "correlated")
oo <- ee$omega_t$cred
ll <- apply(ee$model$lambda, c(3, 4), mean)

expect_equal(as.numeric(oo), c(0.8458204, 0.8838559), tolerance = tol)
expect_equal(ll, matrix(c(0.001341535, -0.001977715, -0.004433084, -0.01164047, 0, 0, -0.02388549,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4025975, 0.5298808,
                          0.4298299, 0.4088522, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.199566,
                          0, 0, 0, 0, 0, 0.4158991, 0.4981815, 0.4208077, 0.5908278, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0.005519423, 0, 0, 0, 0, 0, 0, 0, 0, 0.004320349, -0.0002076825,
                          0.0006509845, -0.002636276, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.02787843,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, -0.000498282, -0.001190331, -0.001030626, 0.001875177),
                        20, 5), tolerance = tol)

pp <- pOmegas(ee)
expect_equal(as.numeric(pp), c(0.748, 1.000), tolerance = tol)



# Frequentist omegas are correct with correlated model, missing listwise, crossloadings
mod <- "
f1 =~ U17_r + U22_r + U29_r + U34_r + U19
f2 =~ U4 + U14 + U19 + U27
f3 =~ U6 + U16 + U28 + U48 + U29_r
f4 =~ U23_r + U31_r + U36_r + U46_r + U34_r
f5 =~ U10_r + U20_r + U35_r + U52_r + U19
"
data(upps, package = "Bayesrel")
set.seed(1234)
ee <- omegasCFA(upps, n.factors = 5, model.type = "correlated", missing = "listwise", model = mod)

expect_equal(c(ee$omega_t$est, ee$omega_t$conf),
             c(0.8675540, 0.8495873, 0.8855208), tolerance = tol)
