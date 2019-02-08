z <- rnorm(100)
z <- z/sum(z)
pvals <- pnorm(abs(z),lower.tail = FALSE)*2
var <- 0.1
maf <- rep(0.05, length(z))
tabledata <- readRDS("./tests/testdata.RDS")

context("Check output of key functions")
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

test_that("prediction functions output value between 0 and 1", {

  x <- pred_logit(tabledata, size = mean(tabledata$claimed.cov))
  expect_true(x >= 0 & x <= 1)

  x1 <- pred_na(tabledata)
  expect_true(x1 >= 0 & x1 <= 1)
})


