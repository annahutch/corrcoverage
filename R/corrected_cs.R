#' corrected_cs
#'
#' @param z Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma Correlation matrix of SNPs
#' @param lower Lower threshold
#' @param upper Upper threshold
#' @param desired.cov The desired coverage of the causal variant in the credible set
#' @param acc Accuracy of corrected coverage to desired coverage (default = 0.0005)
#' @param max.iter Maximum iterations (default = 20)
#'
#' @return list of variants in credible set, required threshold, the corrected coverage and the size of the credible set
#' @export
corrected_cs <- function(z, f, N0, N1, Sigma, lower, upper, desired.cov, acc = 0.0005, max.iter = 20){
  # lower <- 2*desired.cov - 1
  # upper <- min(1,desired.cov + 0.05)
  s = N1/(N0+N1) # proportion of cases
  V = 1/(2 * (N0+N1) * f * (1 - f) * s * (1 - s))
  W = 0.2
  r <- W^2/(W^2 + V)
  pp <- ppfunc(z, V = V) # pp of system in question
  muhat = est_mu(z, f, N0, N1)
  nsnps = length(pp)
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i,]) # nsnp zj vectors for each snp considered causal
  # simulate ERR matrix

  ERR = mvtnorm:::rmvnorm(1000,rep(0,ncol(Sigma)),Sigma)
  pp_ERR = function(Zj, nrep = 1000, Sigma){
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, 1000, length(Zj), byrow=TRUE) # matrix of Zj replicated in each row
    zstar = mexp.zm+ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = coloc:::logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp / rowSums(pp.tmp)
  }
  # simulate pp systems
  pps <- mapply(pp_ERR, zj, MoreArgs = list(nrep = 1000, Sigma = Sigma), SIMPLIFY = FALSE)

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
#' @param Sigma Correlation matrix of SNPs
#' @param lower Lower threshold
#' @param upper Upper threshold
#' @param desired.cov The desired coverage of the causal variant in the credible set
#' @param acc Accuracy of corrected coverage to desired coverage (default = 0.0005)
#' @param max.iter Maximum iterations (default = 20)
#'
#' @return list of variants in credible set, required threshold, the corrected coverage and the size of the credible set
#' @export
corrected_cs_bhat <- function(bhat, V, N0, N1, Sigma, lower, upper, desired.cov, acc = 0.0005, max.iter = 20){
  # lower <- 2*desired.cov - 1
  # upper <- min(1,desired.cov + 0.05)
  z = bhat/sqrt(V)
  W = 0.2
  r <- W^2/(W^2 + V)
  pp = ppfunc(z, V = V) # pp of system in question
  muhat = est_mu_bhat(bhat, V, N0, N1, W = 0.2)
  nsnps = length(pp)
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i,]) # nsnp zj vectors for each snp considered causal
  # simulate ERR matrix

  ERR = mvtnorm:::rmvnorm(1000,rep(0,ncol(Sigma)),Sigma)
  pp_ERR = function(Zj, nrep = 1000, Sigma){
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, 1000, length(Zj), byrow=TRUE) # matrix of Zj replicated in each row
    zstar = mexp.zm+ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = coloc:::logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp / rowSums(pp.tmp)
  }
  # simulate pp systems
  pps <- mapply(pp_ERR, zj, MoreArgs = list(nrep = 1000, Sigma = Sigma), SIMPLIFY = FALSE)

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
