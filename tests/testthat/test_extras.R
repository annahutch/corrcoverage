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
  r <- 200 # no. samples
  c <- 100 # no. SNPs
  h <- matrix(sample(0:1, r*c, replace=TRUE), r, c) # phased haplotype matrix
  cor_matrix <- cor2(h)
  expect_true(dim(cor_matrix)[1]==c)
  expect_true(dim(cor_matrix)[2]==c)
  expect_true(all(dplyr::between(cor_matrix,-1,1.0001)))
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
  V <- coloc:::Var.data.cc(maf, N, 0.5)
  expect_true(all(V>=0))
})


