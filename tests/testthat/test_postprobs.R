z <- rnorm(100)
z <- z/sum(z)
var <- 0.1

test_that("post probs sum to 1", {

  pp <- ppfunc(z, V = var)

  expect_that( length(pp), equals(length(z)) )
  expect_that( sum(pp), equals(1) )
})
