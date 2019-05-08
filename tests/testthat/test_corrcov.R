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
  expect_equal(dim(credset(pp, thr = 0.9))[2], 2)
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
  corr <- corrected_cov(mu = mu, V = V, Sigma = sigma, pp0 = pp, thr = thr)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov function returns a probability", {
  corr <- corrcov(z = z, f = maf, N0 = N, N1 = N, Sigma = sigma, thr = thr)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov_bhat function returns a probability", {
  se <- 0.2 # assume all beta hats have same standard error
  bhats <- z*se
  corr <- corrcov_bhat(bhats, V = V, N0 = N, N1 = N, Sigma = sigma, thr = thr)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov_nvar function returns a probability", {
  corr <- corrcov_nvar(z = z, f = maf, N0 = N, N1 = N, Sigma = sigma, nvar = 2, thr = 0.95, W = 0.2, nrep = 1000)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov_nvar_bhat function returns a probability", {
  se <- 0.2 # assume all beta hats have same standard error
  bhats <- z*se
  corr <- corrcov_nvar_bhat(bhats, V = V, N0 = N, N1 = N, Sigma = sigma, nvar = 2, thr = 0.95, W = 0.2, nrep = 1000)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov_CI returns an appropriate confidence interval", {
  skip("takes too long")
  CI <- corrcov_CI(z = z, f = maf, N0 = N, N1 = N, Sigma = sigma)
  expect_true(dplyr::between(CI[[1]],0,1))
  expect_true(dplyr::between(CI[[2]],0,1))
})

test_that("corrcov_CI_bhat returns an appropriate confidence interval", {
  skip("takes too long")
  se <- 0.2 # assume all beta hats have same standard error
  bhats <- z*se
  CI <- corrcov_CI_bhat(bhat = bhats, V = V, N0 = N, N1 = N, Sigma = sigma)
  expect_true(dplyr::between(CI[[1]],0,1))
  expect_true(dplyr::between(CI[[2]],0,1))
})







