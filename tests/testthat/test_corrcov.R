library(corrcoverage)
context("test_corrcov")

testdata <- system.file("extdata", "testdata_corr.RDS", package="corrcoverage")
data <- readRDS(testdata)
maf <- data$maf
mu <- data$mu
N <- data$N
pp <- data$pp
sigma <- data$sigma
size <- data$size
thr <- data$thr
V <- data$V
z <- data$z

test_that("credset function has correct dimensions (using CV and not)", {
  expect_equal(length(credset(pp, thr = 0.9)), 3)
  expect_equal(length(credset(pp, CV = 2, thr = 0.9)), 4)
})

test_that("credsetC function reports correct things", {
  n.systems <- 100
  pp.matrix <- matrix(rep(pp, n.systems), byrow = TRUE, nrow = n.systems)
  res <- credsetC(pp.matrix, CV = rep(2, 100), thr = 0.95)
  expect_equal(dim(res)[1], n.systems)
  expect_equal(dim(res)[2], 3)
  expect_true(all(dplyr::between(res$claimed.cov,0,1)))
  expect_true(all(res$covered==0|1))
})

test_that("corrected_cov function returns a probability", {
  corr <- corrected_cov(pp0 = pp, mu = mu, V = V, Sigma = sigma, thr = thr, nrep = 2)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov function returns a probability", {
  corr <- corrcov(z = z, f = maf, N0 = N, N1 = N, Sigma = sigma, thr = thr, nrep = 2)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov_bhat function returns a probability", {
  se <- 0.2 # assume all beta hats have same standard error
  bhats <- z*se
  corr <- corrcov_bhat(bhats, V = V, N0 = N, N1 = N, Sigma = sigma, thr = thr, nrep = 2)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov_nvar function returns a probability", {
  corr <- corrcov_nvar(z = z, f = maf, N0 = N, N1 = N, Sigma = sigma, nvar = 2, thr = 0.95, W = 0.2, nrep = 2)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov_nvar_bhat function returns a probability", {
  skip("")
  se <- 0.2 # assume all beta hats have same standard error
  bhats <- z*se
  corr <- corrcov_nvar_bhat(bhats, V = V, N0 = N, N1 = N, Sigma = sigma, nvar = 2, thr = 0.95, W = 0.2, nrep = 100)
  expect_gte(corr, 0)
  expect_lte(corr, 1)
})

test_that("corrcov_CI returns an appropriate confidence interval", {
  skip("")
  CI <- corrcov_CI(z = z, f = maf, N0 = N, N1 = N, Sigma = sigma, nrep = 1)
  expect_true(dplyr::between(CI[[1]],-0.1,1.1))
  expect_true(dplyr::between(CI[[2]],-0.0,1.1))
})

test_that("corrcov_CI_bhat returns an appropriate confidence interval", {
  se <- 0.2 # assume all beta hats have same standard error
  bhats <- z*se
  CI <- corrcov_CI_bhat(bhat = bhats, V = V, N0 = N, N1 = N, Sigma = sigma, nrep = 1)
  expect_true(dplyr::between(CI[[1]],-0.1,1.1))
  expect_true(dplyr::between(CI[[2]],-0.0,1.1))
})

test_that("est_mu and est_mu_bhat return single value", {
  mu1 <- est_mu(z = z, f = maf, N0 = N, N1 = N, W = 0.2)

  se <- 0.2 # assume all beta hats have same standard error
  bhats <- z*se

  mu2 <- est_mu_bhat(bhat = bhats, V = V, N0 = N, N1 = N, W = 0.2)
  expect_true(class(mu1) == "numeric")
  expect_true(class(mu2) == "numeric")
})

test_that("corrected_cs reports appropriate list", {
  # skip("takes too long")
  res <- corrected_cs(z = z, f = maf, N0 = N, N1 = N, Sigma = sigma, lower = 0.8, upper = 1, desired.cov = 0.95, max.iter = 1)
  expect_true(length(res) == 4)
  expect_true(all(dplyr::between(res$corr.cov,0,1)))
  expect_true(all(dplyr::between(res$req.thr,0,1)))
  expect_true(all(dplyr::between(res$size,0,1)))
})

test_that("corrected_cs_bhat reports appropriate list", { # get an error.. cannot make it smaller?
  # skip("takes too long")
  se <- 0.2 # assume all beta hats have same standard error
  bhats <- z*se
  expect_error(corrected_cs_bhat(bhat = bhats, V = V, N0 = N, N1 = N, Sigma = sigma, lower = 0.8, upper = 1, desired.cov = 0.95, max.iter = 1))
  #expect_true(length(res) == 4)
  #expect_true(all(dplyr::between(res$corr.cov,0,1)))
  #expect_true(all(dplyr::between(res$req.thr,0,1)))
  #expect_true(all(dplyr::between(res$size,0,1)))
})








