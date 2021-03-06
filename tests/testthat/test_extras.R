library(corrcoverage)
context("test_extras")

N <- 5000
z <- rnorm(100) # z-scores
pvals <- pnorm(abs(z),lower.tail = FALSE)*2 # p-values
var <- rep(0.01, length(z)) # variance of estimated effect size
maf <- rep(0.05, length(z)) # minor allele freqs
r <- 200 # no. samples
c <- 100 # no. SNPs
h <- matrix(sample(0:1, r*c, replace=TRUE), r, c) # phased haplotype matrix

testdata <- system.file("extdata", "testdata.RDS", package="corrcoverage")
data <- readRDS(testdata)

test_that("cor2 finds correlation matrix", {
  r <- 10 # no. samples
  c <- 100 # no. SNPs
  h <- matrix(sample(0:1, r*c, replace=TRUE), r, c) # phased haplotype matrix
  cor_matrix <- cor2(h)
  expect_true(dim(cor_matrix)[1]==c)
  expect_true(dim(cor_matrix)[2]==c)
})

test_that("z_sim simulates the correct number of marginal Z scores", {
  zj = rep(0, 100)
  zj[2] = 5
  res = z_sim(Zj = zj, Sigma = cor2(h), nrep = 10)
  expect_true(dim(res)[1] == 10)
  expect_true(dim(res)[2] == 100)
})

test_that("zj_pp returns probabilities and has correct dimensions", {
  zj = rep(0, 100)
  zj[2] = 5
  res = zj_pp(Zj = zj, V = var, nrep = 10, W = 0.2, Sigma = cor2(h))
  expect_true(all( res > 0 & res < 1))
  expect_true(dim(res)[1] == 10)
  expect_true(dim(res)[2] == 100)
})

test_that("length of posterior probs is the same as input (or plus 1 if null model considered)", {

  pvals.pp.null <- pvals_pp(pvals, f = maf, type = "cc", N = 10000, s = 0.5)
  expect_that( length(pvals.pp.null), equals(length(z)+1))
  expect_true(all(pvals.pp.null>=0))

  pp <- ppfunc(z, V = var)
  expect_that( length(pp), equals(length(z)))
  expect_true(all(pp>=0))

  z.pp.null <- z0_pp(z, f = maf, type = "cc", N = 10000, s = 0.5)
  expect_that( length(z.pp.null), equals(length(z)+1))
  expect_true(all(z.pp.null>=0))

})

test_that("no missing values", {
  pvals.pp.null <- pvals_pp(pvals, f = maf, type = "cc", N = 10000, s = 0.5)
  expect_identical(pvals.pp.null, na.omit(pvals.pp.null))

  pp <- ppfunc(z, V = var)
  expect_identical(pp, na.omit(pp))

  z.pp.null <- z0_pp(z, f = maf, type = "cc", N = 10000, s = 0.5)
  expect_identical(z.pp.null, na.omit(z.pp.null))
})

test_that("check V in ppfunc", {
  expect_error(ppfunc(z))
  V <- Var.data.cc(f = maf, N, 0.5)
  expect_true(all(V>=0))
})

test_that("ppfunc.mat returns probabilities", {
  V <- Var.data.cc(maf, N, 0.5)
  zstar <- matrix(rnorm(10000), nrow = 100)
  res <- ppfunc.mat(zstar, V = V)
  expect_equal(sum(res[1,]), 1)
  expect_true(all(dplyr::between(res,0,1)))
})

test_that("approx.bf.p returns 4 columns", {
  res <- approx.bf.p(pvals, f = maf, type = "cc", N = N, s = 0.5, W = 0.2)
  expect_true(dim(res)[2]==4)
})

test_that("bf_func returns vector of correct length", {
  varbeta <- Var.data.cc(maf, N, 0.5)
  expect_true(length(bf_func(z = z, V = varbeta))==length(z))
})

test_that("calculating bfs over a matrix performs row wise", {
  W = 0.2
  V = Var.data.cc(maf, N, 0.5)
  r = W^2/(W^2 + V)
  zj = rep(0, 100)
  zj[2] = 5
  nrep = 10
  Z = z_sim(Zj = zj, Sigma = cor2(h), nrep = nrep)
  bf_matrix = 0.5 * t(log(1 - r) + (r * t(Z^2)))
  bf_apply = sapply(1:nrep, function(i) bf = 0.5 * (log(1 - r) + (r * Z[i,]^2))) %>% t()
  expect_identical(bf_matrix, bf_apply)
})

test_that("bf_func only accepts parameters of the same class", {
  V <- Var.data.cc(f = maf, N, 0.5)
  V_matrix = t(replicate(2, V))
  z_matrix = t(replicate(2, z))
  expect_error(bf_func(z = z, V = V_matrix))
  expect_error(bf_func(z = z_matrix, V = V))
})

test_that(".zj_abf works", {
  W = 0.2
  V = Var.data.cc(f = maf, N, 0.5)
  r = W^2/(W^2+V)
  nrep = 10
  sigma =  cor2(h)
  ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(sigma)), sigma)
  res <- .zj_abf(Zj = z, sigma, nrep, ERR, r)
  expect_identical(dim(res), dim(ERR))
})

test_that(".zj_abf only accepts parameters of correct class", {
  W = 0.2
  V = Var.data.cc(f = maf, N, 0.5)
  r = W^2/(W^2+V)
  nrep = 10
  sigma =  cor2(h)
  z_wrong = t(replicate(2, z))
  r_wrong = t(replicate(2, r))
  ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(sigma)), sigma)
  expect_error(.zj_abf(Zj = z_wrong, sigma, nrep, ERR, r))
  expect_error(.zj_abf(Zj = z, sigma, nrep, ERR[1,], r))
  expect_error(.zj_abf(Zj = z, sigma, nrep, ERR, r_wrong))
})

test_that(".zj_pp only accepts parameters of correct class", {
  W = 0.2
  V = Var.data.cc(f = maf, N, 0.5)
  r = W^2/(W^2+V)
  nrep = 10
  sigma =  cor2(h)
  z_wrong = t(replicate(2, z))
  r_wrong = t(replicate(2, r))
  ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(sigma)), sigma)
  expect_error(.zj_pp(Zj = z_wrong, sigma, nrep, ERR, r))
  expect_error(.zj_pp(Zj = z, sigma, nrep, ERR[1,], r))
  expect_error(.zj_pp(Zj = z, sigma, nrep, ERR, r_wrong))
})
