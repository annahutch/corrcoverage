#' Variance of the estimated effect size for case-control data
#'
#' @title Variance of the estimated effect size for case-control data
#' @param f Minor allele frequencies
#' @param N Total sample size (N0+N1)
#' @param s Proportion of cases (N1/N0+N1)
#' @return Variance of estimated effect size \eqn{\hat{\beta}}, V.
#' @author Claudia Giambartolomei
#'
#' @examples
#'
#' maf =  runif(100, 0.05, 0.5)
#' N0 = 5000 # number of controls
#' N1 = 5000 # number of cases
#'
#' Var.data.cc(f = maf, N = N0 + N1, s = N1/(N0+N1))
#'
#' @export
Var.data.cc <- function(f, N, s) {
    1/(2 * N * f * (1 - f) * s * (1 - s))
}

##' Internal function, logsum
##'
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##' @title logsum
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
##' @author Claudia Giambartolomei
##' @export
logsum <- function(x) {
    my.max <- max(x)
    my.res <- my.max + log(sum(exp(x - my.max)))
    return(my.res)
}

#' Correlation matrix of SNPs
#'
#' Quick function to find a correlation matrix
#' @title Correlation matrix of SNPS
#' @param x Phased haplotype matrix, rows as samples and columns as SNPs
#' @return Correlation matrix
#' @author Chris Wallace
#' @export
cor2 <- function(x) {
    1/(NROW(x) - 1) * crossprod(scale(x, TRUE, TRUE))
}

#' Proportion of simulated credible sets containing the causal variant
#'
#' @rdname prop_cov
#' @title Proportion of credible sets containing the causal variant
#' @param x data.frame with a binary 'covered' column
#' @return Proportion of x with x = 1
#' @author Anna Hutchinson
prop_cov <- function(x) {
    mean(x$covered)
}

#' @title Estimate the true effect at the causal variant using Z-scores and MAFs
#'
#' @param z Vector of marginal Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param W Prior for the standard deviation of the effect size parameter, beta, default 0.2
#'
#' @return Estimate of the true effect at the causal variant
#'
#' @examples
#'
#' nsnps <- 100
#' z_scores <- rnorm(nsnps, 0, 3) # simulate a vector of Z-scores
#' N0 <- 5000 # number of controls
#' N1 <- 5000 # number of cases
#'
#' maf <- runif(nsnps, 0.05, 0.5)
#'
#' est_mu(z = z_scores, f = maf, N0 = N0, N1 = N1)
#'
#' @export
#'
#' @author Anna Hutchinson
est_mu <- function(z, f, N0, N1, W = 0.2) {
  stopifnot(class(z)=="numeric") # ensure z is not a matrix of simulated z scores
  V = 1/(2 * (N0 + N1) * f * (1 - f) * (N1/(N0 + N1)) * (1 - (N1/(N0 + N1))))
  r = W^2/(W^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  denom = logsum(lABF)
  pp = exp(lABF - denom)  # convert back from log scale
  sum(abs(z)*pp)
}

#' @title Estimate the true effect at the causal variant using estimated
#' effect sizes and their standard errors
#'
#' @param bhat Vector of estimated effect sizes
#' @param V Prior variance for estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param p1 Prior probability a SNP is associated with the trait, default 1e-4
#' @param W Prior for the standard deviation of the effect size parameter, beta
#'
#' @return Estimate of the true effect at the causal variant
#'
#' @examples
#'
#' nsnps <- 100
#' N0 <- 5000 # number of controls
#' N1 <- 5000 # number of cases
#'
#' maf <- runif(nsnps, 0.05, 0.3)
#'
#' varbeta <- Var.data.cc(f = maf, N = N0 + N1, s = N1/(N0+N1))
#'
#' bhats = rnorm(nsnps,0,0.2) # log(OR)
#'
#' est_mu_bhat(bhat = bhats, V = varbeta, N0 = N0, N1 = N1)
#'
#' @export
#'
#' @author Anna Hutchinson
est_mu_bhat <- function(bhat, V, N0, N1, p1 = 1e-4, W = 0.2) {
    stopifnot(class(bhat)=="numeric") # ensure bhat is not a matrix of simulated z scores
    z = bhat/sqrt(V)
    r = W^2/(W^2 + V)
    lABF = 0.5 * (log(1 - r) + (r * z^2))
    denom = logsum(lABF)
    pp = exp(lABF - denom)  # convert back from log scale
    sum(abs(z)*pp)
}

#' Credible set of putative causal variants
#'
#' If the CV parameter is supplied (index of causal variant) then the
#' output includes a binary indicator of whether the CV is contained in the set
#' @title Credible set of genetic variants
#' @param pp Vector of posterior probabilities of causality
#' @param CV Optional parameter: Index of CV
#' @param thr Minimum threshold for credible set size (default is 0.95)
#' @return list of the variants in the credible set, the claimed.cov (cumulative sum of the posterior probabilities of the variants forming the credible set), binary covered indicator (1 if CV is contained in the credible set) and nvar (number of variants in the set)
#'
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#' pp <- rnorm(nsnps, 0.3, 0.05)
#' pp <- pp/sum(pp)
#'
#' credset(pp, thr = 0.9)
#'
#' iCV <- 71
#'
#' credset(pp, CV = iCV, thr = 0.9)
#'
#' @export
#' @author Anna Hutchinson
credset <- function(pp, CV, thr = 0.95) {
  o = order(pp, decreasing = TRUE)  # order index for true pp
  cumpp = cumsum(pp[o])  # cum sums of ordered pps
  wh = which(cumpp > thr)[1]  # how many needed to exceed thr
  names(wh) = NULL
  size = cumpp[wh]
  names(size) = NULL
  if (missing(CV)) {
    list(credset = o[1:wh], claimed.cov = size, nvar = wh)
  } else {
    contained = as.numeric(CV %in% o[1:wh])
    list(credset = o[1:wh], claimed.cov = size, covered = contained, nvar = wh)
  }
}

#' Quick credset function for matrix of posterior probabilities (using RCpp)
#'
#' @title Credible set of variants from matrix of PPs
#' @param pp Matrix of posterior probabilities of causality (one row per system)
#' @param CV Vector of CV indices (one per system/row)
#' @param thr Minimum threshold for credible set size (default is 0.95)
#'
#' @return Data.frame of claimed coverage (sum of posterior probabilities of variants in the set), binary covered indicator and number of variants (nvar).
#' @useDynLib corrcoverage
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#'
#' # simulate matrix of posterior probabilities
#' # 1 simulation per row
#'
#' pp <- matrix(rnorm(nsnps*100, 0.3, 0.05), ncol = nsnps)
#' pp <- pp/rowSums(pp)
#'
#' iCV <- 71
#'
#' # credsetC(pp, CV = iCV, thr = 0.9)
#'
#'
#' @export
credsetC <- function(pp, CV, thr = 0.95) {
  ret = credsetmat(pp, CV, thr)  ## list 1 = wh, 2 = size, 3=contained
  data.frame(claimed.cov = ret[[2]], covered = ret[[3]], nvar = ret[[1]])
}

#' Calculate ABFs from Z scores
#'
#' Note, z and V should both be vectors or both be matrices
#' @title Calculate ABFs from Z scores
#' @param z Vector of Z-scores
#' @param V Variance of the estimated effect size
#' @param W Prior for the standard deviation of the effect size parameter, beta (default 0.2)
#'
#' @return ABFs
#'
#' @examples
#'
#' set.seed(1)
#' nsnps = 100
#' N0 = 5000
#' N1 = 5000
#' z_scores <- rnorm(nsnps, 0, 3)
#'
#' ## generate example LD matrix and MAFs
#' library(mvtnorm)
#' nsamples = 1000
#'
#' simx <- function(nsnps, nsamples, S, maf=0.1) {
#'     mu <- rep(0,nsnps)
#'     rawvars <- rmvnorm(n=nsamples, mean=mu, sigma=S)
#'     pvars <- pnorm(rawvars)
#'     x <- qbinom(1-pvars, 1, maf)
#'}
#'
#' S <- (1 - (abs(outer(1:nsnps,1:nsnps,`-`))/nsnps))^4
#' X <- simx(nsnps,nsamples,S)
#' maf <- colMeans(X)
#'
#' varbeta = Var.data.cc(f = maf, N = N0 + N1, s = 0.5)
#'
#' bf_func(z_scores, V = varbeta)
#'
#' @export
bf_func <- function(z, V, W = 0.2){
  stopifnot(class(z)==class(V))
  r = W^2/(W^2 + V)
  0.5 * (log(1 - r) + (r * z^2))
}
