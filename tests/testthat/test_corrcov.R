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

test_that("credset function", {
  expect_equal(dim(credset(pp, thr = 0.9))[2], 2)
  expect_equal(dim(credset(pp, CV = 2, thr = 0.9))[2], 3)
})

test_that("corrected_cov function returns a probability", {
  skip('skip')
  corr <- corrected_cov(mu = mu, V = V, Sigma = sigma, pp0 = pp, thresh = thr, size = size)
  expect_true(corr>=0 & corr<=1)
})

test_that("corrcov function returns a probability", {
  skip('skip')
  corr <- corrcov(z = z, f = maf, N0 = N, N1 = N, Sigma = sigma, thr = thr)
  expect_true(corr>=0 & corr<=1)
})




