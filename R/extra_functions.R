#' Variance of the estimated effect size for case-control data
#'
#' @title Variance of the estimated effect size for case-control data
#' @param f Minor allele frequencies
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1)
#' @return Variance of estimated effect size \eqn{\hat{\beta}}, V.
#' @author Claudia Giambartolomei
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
logsum <- function(x) {
    my.max <- max(x)
    my.res <- my.max + log(sum(exp(x - my.max)))
    return(my.res)
}

#' Correlation matrix of SNPS
#'
#' A quick function to find a correlation matrix
#' @title Correlation matrix of SNPS
#' @param x Phased haplotype matrix, rows as samples and columns as SNPs
#' @return Correlation matrix
#' @export
#' @author Chris Wallace
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
#' @export
prop_cov <- function(x) {
    sum(x$covered)/length(x$covered)
}

#' @title Estimate the true effect at the causal variant
#'
#' @param z Vector of marginal Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param W Prior for the standard deviation of the effect size parameter beta
#'
#' @return Estimate of the true effect at the causal variant
#' @export
#'
#' @author Anna Hutchinson
est_mu <- function(z, f, N0, N1, W = 0.2) {
    V = 1/(2 * (N0 + N1) * f * (1 - f) * (N1/(N0 + N1)) * (1 - (N1/(N0 + N1))))
    r = W^2/(W^2 + V)
    lABF = 0.5 * (log(1 - r) + (r * z^2))
    p1 = 1e-04  # hard coded
    nsnps = length(lABF)
    prior = c(1 - nsnps * p1, rep(p1, nsnps))
    tmp = c(1, lABF)  # add on extra for null model
    my.denom = logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    ph0.tmp = tmp1/sum(tmp1)

    ph0 = ph0.tmp[1]  # prob of the null
    mean(c(sum(abs(z) * ph0.tmp[-1]), (1 - ph0.tmp[1]) * max(abs(z))))
}

#' @title Estimate the true effect at the causal variant
#'
#' @param bhat Vector of estimated effect sizes
#' @param V Prior variance for estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param W Prior for the standard deviation of the effect size parameter beta
#'
#' @return Estimate of the true effect at the causal variant
#' @export
#'
#' @author Anna Hutchinson
est_mu_bhat <- function(bhat, V, N0, N1, W = 0.2) {
    z = bhat/sqrt(V)
    r = W^2/(W^2 + V)
    lABF = 0.5 * (log(1 - r) + (r * z^2))
    p1 = 1e-04  # hard coded
    nsnps = length(lABF)
    prior = c(1 - nsnps * p1, rep(p1, nsnps))
    tmp = c(1, lABF)  # add on extra for null model
    my.denom = logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    ph0.tmp = tmp1/sum(tmp1)
    ph0 = ph0.tmp[1]  # prob of the null
    mean(c(sum(abs(z) * ph0.tmp[-1]), (1 - ph0.tmp[1]) * max(abs(z))))
}

#' Get credible set of variants
#'
#' If the CV parameter is supplied (index of causal variant) then the
#' output includes a binary indicator of whether the CV is contained in the set
#' @title Get credible set of variants
#' @param pp Vector of posterior probabilities of causality
#' @param CV Optional parameter: Index of CV
#' @param thr Minimum threshold for credible set size (default is 0.95)
#' @export
#' @return list of the variants in the credible set, the claimed.cov (cumulative sum of the posterior probabilities of the variants forming the credible set), binary covered indicator (1 if CV is contained in the credible set) and nvar (number of variants in the set)
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

#' Quicker credset function for matrix of posterior probabilities (using RCpp)
#'
#' @title Get credible set of variants from matrix of pps (Rcpp)
#' @param pp Matrix of posterior probabilities of causality (one row per system)
#' @param CV Vector of CV indices (one per system/row)
#' @param thr Minimum threshold for credible set size (default is 0.95)
#'
#' @return Data.frame of claimed coverage (sum of posterior probabilities of variants in the set), binary covered indicator and number of variants (nvar).
#' @useDynLib corrcoverage
#' @importFrom Rcpp sourceCpp
#' @export
credsetC <- function(pp, CV, thr = 0.95) {
  ret = credsetmat(pp, CV, thr)  ## list 1 = wh, 2 = size, 3=contained
  data.frame(claimed.cov = ret[[2]], covered = ret[[3]], nvar = ret[[1]])
}
