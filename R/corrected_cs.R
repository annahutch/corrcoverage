#' corrected_cs
#'
#' @param z Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param lower Lower threshold
#' @param upper Upper threshold
#' @param desired.cov The desired coverage of the causal variant in the credible set
#' @param acc Accuracy of corrected coverage to desired coverage (default = 0.0005)
#' @param max.iter Maximum iterations (default = 20)
#'
#' @return list of variants in credible set, required threshold, the corrected coverage and the size of the credible set
#' @export
corrected_cs <- function(z, f, N0, N1, lower, upper, desired.cov, acc = 0.0005, max.iter = 20){
  # lower <- 2*desired.cov - 1
  # upper <- min(1,desired.cov + 0.05)
  s = N1/(N0+N1) # proportion of cases
  V = 1/(2 * (N0+N1) * f * (1 - f) * s * (1 - s))
  pp <- ppfunc(z, V = V) # pp of system in question
  muhat = est_mu(z, f, N0, N1)
  nsnps = length(pp)
  temp <- diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj <- do.call(c, apply(temp, 1, list))  # nsnp zj vectors for each snp considered causal

  # simulate pp systems
  pps <- mapply(zj_pp, zj, V = V, MoreArgs = list(nrep = 1000, Sigma = LD), SIMPLIFY = FALSE)

  n_pps <- length(pps)
  args <- 1:length(pp)

  f <- function(thr){ # finds the difference between corrcov and desired.cov
    d5 <- lapply(1:n_pps, function(x) {
      credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })
    prop_cov <- lapply(d5, prop_cov) %>% unlist()
    sum(prop_cov * pp) - desired.cov
  }

  # initalize
  N=1
  fa = f(lower)
  fb = f(upper)

  if (fa * fb > 0) {
    stop("No root in range, increase window")
  } else {
    fc = min(fa, fb)
    while (N < max.iter & !dplyr::between(fc, 0, acc)) {
      c = lower + (upper-lower)/2
      fc = f(c)
      print(paste("thr: ", c, ", cov: ", desired.cov + fc))

      if (fa * fc < 0) {
        upper = c
        fb = fc
      } else if (f(upper) * fc < 0) {
        lower = c
        fa = fc
      }
      N = N + 1
    }
    # df <- data.frame(req.thr = c, corr.cov = desired.cov + fc)
  }
  o <- order(pp, decreasing = TRUE)  # order index for true pp
  cumpp <- cumsum(pp[o])  # cum sums of ordered pps
  wh <- which(cumpp > c)[1]  # how many needed to exceed thr
  list(credset = names(pp)[o[1:wh]], req.thr = c, corr.cov = desired.cov + fc, size = cumpp[wh])
}

#' corrected_cs_bhat
#'
#' @param bhat Estimated effect sizes
#' @param V Prior variance of estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param lower Lower threshold
#' @param upper Upper threshold
#' @param desired.cov The desired coverage of the causal variant in the credible set
#' @param acc Accuracy of corrected coverage to desired coverage (default = 0.0005)
#' @param max.iter Maximum iterations (default = 20)
#'
#' @return list of variants in credible set, required threshold, the corrected coverage and the size of the credible set
#' @export
corrected_cs_bhat <- function(bhat, V, N0, N1, lower, upper, desired.cov, acc = 0.0005, max.iter = 20){
  # lower <- 2*desired.cov - 1
  # upper <- min(1,desired.cov + 0.05)
  z = bhat/sqrt(V)
  pp = ppfunc(z, V = V) # pp of system in question
  muhat = est_mu(z, f, N0, N1)
  nsnps = length(pp)
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = do.call(c, apply(temp, 1, list))  # nsnp zj vectors for each snp considered causal

  # simulate pp systems
  pps <- mapply(zj_pp, zj, V = V, MoreArgs = list(nrep = 1000, Sigma = LD), SIMPLIFY = FALSE)

  n_pps <- length(pps)
  args <- 1:length(pp)

  f <- function(thr){ # finds the difference between corrcov and desired.cov
    d5 <- lapply(1:n_pps, function(x) {
      credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })
    prop_cov <- lapply(d5, prop_cov) %>% unlist()
    sum(prop_cov * pp) - desired.cov
  }

  # initalize
  N=1
  fa = f(lower)
  fb = f(upper)

  if (fa * fb > 0) {
    stop("No root in range, increase window")
  } else {
    fc = min(fa, fb)
    while (N < max.iter & !dplyr::between(fc, 0, acc)) {
      c = lower + (upper-lower)/2
      fc = f(c)
      print(paste("thr: ", c, ", cov: ", desired.cov + fc))

      if (fa * fc < 0) {
        upper = c
        fb = fc
      } else if (f(upper) * fc < 0) {
        lower = c
        fa = fc
      }
      N = N + 1
    }
    # df <- data.frame(req.thr = c, corr.cov = desired.cov + fc)
  }
  o <- order(pp, decreasing = TRUE)  # order index for true pp
  cumpp <- cumsum(pp[o])  # cum sums of ordered pps
  wh <- which(cumpp > c)[1]  # how many needed to exceed thr
  list(credset = names(pp)[o[1:wh]], req.thr = c, corr.cov = desired.cov + fc, size = cumpp[wh])
}
