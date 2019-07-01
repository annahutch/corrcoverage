#' @title Corrected credible set with desired coverage of the CV
#'
#' @param z Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma Correlation matrix of SNPs
#' @param lower Lower threshold (default = 0)
#' @param upper Upper threshold (default = 1)
#' @param desired.cov The desired coverage of the causal variant in the credible set
#' @param acc Accuracy of corrected coverage to desired coverage (default = 0.005)
#' @param max.iter Maximum iterations (default = 20)
#' @return List of variants in credible set, required threshold, the corrected coverage and the size of the credible set
#'
#' @examples
#'
#' # In this example, the function is used to find a corrected 95% credible set
#' # using Z-scores and MAFs, that is the smallest set of variants
#' # required such that the resultant credible set has coverage close to (/within
#' # some accuracy of) the "desired coverage" (here set to 0.95). Max.iter parameter
#' # defines the maximum number of iterations to try in the root bisection algorithm,
#' # this should be increased to ensure convergence to the desired coverage, but is set
#' # to 2 here for speed (and thus the resultant credible set will not be accurate).
#'
#' set.seed(2)
#' nsnps = 100
#' N0 = 5000
#' N1 = 5000
#' z_scores <- rnorm(nsnps, 0, 3) # simulate a vector of Z-scores
#'
#' # simulate fake haplotypes to obtain MAFs and LD matrix
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' maf <- colMeans(haps)
#' LD <- cor2(haps)
#'
#' corrected_cs(z = z_scores, f = maf, N0, N1, Sigma = LD, desired.cov = 0.95, max.iter = 2)
#' # max.iter set low for speed, should be set to at least
#' # the default to ensure convergence to desired coverage
#'
#' @export
#' @author Anna Hutchinson
corrected_cs <- function(z, f, N0, N1, Sigma, lower = 0, upper = 1, desired.cov, acc = 0.005, max.iter = 20){
  s = N1/(N0+N1) # proportion of cases
  V = 1/(2 * (N0+N1) * f * (1 - f) * s * (1 - s))
  W = 0.2
  r = W^2/(W^2 + V)
  pp = ppfunc(z, V = V) # pp of system in question
  muhat = est_mu(z, f, N0, N1)
  nsnps = length(pp)
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i,]) # nsnp zj vectors for each snp considered causal
  nrep = 1000

  # simulate ERR matrix
  ERR = mvtnorm::rmvnorm(nrep, rep(0,ncol(Sigma)), Sigma)
  pp_ERR = function(Zj){
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow=TRUE) # matrix of Zj replicated in each row
    zstar = mexp.zm+ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = logsum(bf)
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp / rowSums(pp.tmp)
  }
  # simulate pp systems
  pps = mapply(pp_ERR, zj, SIMPLIFY = FALSE)

  n_pps = length(pps)
  args = 1:length(pp)

  f <- function(thr){ # finds the difference between corrcov and desired.cov
    d5 <- lapply(1:n_pps, function(x) {
      credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })
    prop_cov <- lapply(d5, prop_cov) %>% unlist()
    sum(prop_cov * pp) - desired.cov
  }

  o = order(pp, decreasing = TRUE)  # order index for true pp
  cumpp = cumsum(pp[o])  # cum sums of ordered pps

  corrcov.tmp <- f(desired.cov)
  nvar.tmp <- which(cumpp > desired.cov)[1]

  if(corrcov.tmp > 0 & nvar.tmp == 1) stop("Cannot make credible set smaller")

  # initalize
  N=1
  fa = f(lower)
  fb = f(upper)

  if (fa * fb > 0) {
    stop("No root in range, increase window")
  } else {
    fc = min(fa, fb)
    while (N <= max.iter & !dplyr::between(fc, 0, acc)) {
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
  }
  wh = which(cumpp > c)[1]  # how many needed to exceed thr
  size = cumpp[wh]
  names(size) = NULL
  list(credset = names(pp)[o[1:wh]], req.thr = c, corr.cov = desired.cov + fc, size = size)
}

#' @title Get new credible set with desired coverage of the CV
#'
#' @param bhat Estimated effect sizes
#' @param V Prior variance of estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma Correlation matrix of SNPs
#' @param lower Lower threshold (default = 0)
#' @param upper Upper threshold (default = 1)
#' @param desired.cov The desired coverage of the causal variant in the credible set
#' @param acc Accuracy of corrected coverage to desired coverage (default = 0.005)
#' @param max.iter Maximum iterations (default = 20)
#'
#' @return List of variants in credible set, required threshold, the corrected coverage and the size of the credible set
#'
#' @examples
#'
#' # In this example, the function is used to find a corrected 95% credible set
#' # using bhats and their standard errors, that is the smallest set of variants
#' # required such that the resultant credible set has coverage close to (/within
#' # some accuracy of) the "desired coverage" (here set to 0.95). Max.iter parameter
#' # defines the maximum number of iterations to try in the root bisection algorithm,
#' # this should be increased to ensure convergence to the desired coverage, but is set
#' # to 2 here for speed (and thus the resultant credible set will not be accurate).
#'
#' set.seed(1)
#' nsnps <- 100
#' iCV <- 4
#' N0 <- 5000 # number of controls
#' N1 <- 5000 # number of cases
#'
#' # simulate fake haplotypes to obtain MAFs and LD matrix
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' maf <- colMeans(haps)
#' LD <- cor2(haps)
#'
#' varbeta <- Var.data.cc(f = maf, N = N0 + N1, s = N1/(N0+N1))
#'
#' beta <- rep(0, nsnps)
#' beta[iCV] <- 5
#'
#' bhats <- rnorm(beta, varbeta)
#'
#' corrcov_cs_bhat(bhat = bhats, V = varbeta, N0, N1, Sigma = LD, desired.cov = 0.95, max.iter = 2)
#' # max.iter set low for speed, should be set to at
#' # least the default to ensure convergence to desired coverage
#'
#' @export
#' @author Anna Hutchinson
corrected_cs_bhat <- function(bhat, V, N0, N1, Sigma, lower = 0, upper = 1, desired.cov, acc = 0.005, max.iter = 20){
  z = bhat/sqrt(V)
  W = 0.2
  r = W^2/(W^2 + V)
  pp = ppfunc(z, V = V) # pp of system in question
  muhat = est_mu_bhat(bhat, V, N0, N1, W = 0.2)
  nsnps = length(pp)
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i,]) # nsnp zj vectors for each snp considered causal
  nrep = 1000

  # simulate ERR matrix
  ERR = mvtnorm::rmvnorm(nrep, rep(0,ncol(Sigma)), Sigma)
  pp_ERR = function(Zj){
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow=TRUE) # matrix of Zj replicated in each row
    zstar = mexp.zm+ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = logsum(bf)
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp / rowSums(pp.tmp)
  }
  # simulate pp systems
  pps = mapply(pp_ERR, zj, SIMPLIFY = FALSE)

  n_pps = length(pps)
  args = 1:length(pp)

  f <- function(thr){ # finds the difference between corrcov and desired.cov
    d5 <- lapply(1:n_pps, function(x) {
      credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })
    prop_cov <- lapply(d5, prop_cov) %>% unlist()
    sum(prop_cov * pp) - desired.cov
  }

  o = order(pp, decreasing = TRUE)  # order index for true pp
  cumpp = cumsum(pp[o])  # cum sums of ordered pps

  corrcov.tmp <- f(desired.cov)
  nvar.tmp <- which(cumpp > desired.cov)[1]

  if(corrcov.tmp > 0 & nvar.tmp == 1) stop("Cannot make credible set smaller")

  # initalize
  N=1
  fa = f(lower)
  fb = f(upper)

  if (fa * fb > 0) {
    stop("No root in range, increase window")
  } else {
    fc = min(fa, fb)
    while (N <= max.iter & !dplyr::between(fc, 0, acc)) {
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
  }
  wh = which(cumpp > c)[1]  # how many needed to exceed thr

  size = cumpp[wh]
  names(size) = NULL
  list(credset = names(pp)[o[1:wh]], req.thr = c, corr.cov = desired.cov + fc, size = size)
}
