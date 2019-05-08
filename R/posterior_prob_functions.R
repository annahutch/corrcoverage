#' Wakefield's log asymptotic Bayes factor (lABF) with prior standard deviation of effect size as a parameter
#'
#' ([Wakefield et al. 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359)
#' This function converts p-values to log ABFs, also reporting intermediate calculations
#' @title Find approx. Bayes factors
#' @param pvals p-values
#' @param f Minor allele frequencies
#' @param type Type of experiment ('quant' or 'cc')
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=='quant'
#' @param W Prior for the standard deviation of the effect size parameter beta (W=0.2 default)
#' @return data.frame containing lABF and intermediate calculations
#' @export
approx.bf.p <- function(pvals, f, type, N, s, W = 0.2) {
    V = 1/(2 * N * f * (1 - f) * s * (1 - s))
    z = stats::qnorm(0.5 * pvals, lower.tail = FALSE)
    r = W^2/(W^2 + V)
    lABF = 0.5 * (log(1 - r) + (r * z^2))
    ret = data.frame(V, z, r, lABF)
    return(ret)
}

#' Posterior probabilities of causality from p-values
#'
#' This function converts p-values to posterior probabilities of causality, including the null model of no genetic effects
#' @title Find posterior probabilities of causality
#' @param pvals p-values of SNPs
#' @param f Minor allele frequencies
#' @param type Type of experiment ('quant' or 'cc')
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=='quant'
#' @param W Prior for the standard deviation of the effect size parameter, beta (W = 0.2 default)
#' @return Posterior probabilities of null model (no genetic effect) and causality for each SNP
#' @export
#' @author Anna Hutchinson
pvals_pp <- function(pvals, f, type, N, s, W = 0.2) {
    V = 1/(2 * N * f * (1 - f) * s * (1 - s))
    z = stats::qnorm(0.5 * pvals, lower.tail = FALSE)
    r = W^2/(W^2 + V)
    lABF = 0.5 * (log(1 - r) + (r * z^2))
    p1 = 1e-04  # hard coded
    nsnps = length(lABF)
    prior = c(1 - nsnps * p1, rep(p1, nsnps))
    tmp = c(1, lABF)  # add on extra for null model
    my.denom = coloc:::logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    tmp1/sum(tmp1)
}

#' Posterior probabilities of causality from marginal Z-scores with prior SD as a parameter
#'
#' This function converts Z-scores to posterior probabilities of causality, including the null model of no genetic effects
#' @title Find posterior probabilities of causality
#' @param z Marginal Z-scores of SNPs
#' @param f Minor allele frequencies
#' @param type Type of experiment ('quant' or 'cc')
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=='quant'
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @return Posterior probabilities of null model (no genetic effect) and causality for each SNP
#' @export
#' @author Anna Hutchinson
z0_pp <- function(z, f, type, N, s, W = 0.2) {
    V = 1/(2 * N * f * (1 - f) * s * (1 - s))
    r = W^2/(W^2 + V)
    lABF = 0.5 * (log(1 - r) + (r * z^2))
    p1 = 1e-04  # hard coded
    nsnps = length(lABF)
    prior = c(1 - nsnps * p1, rep(p1, nsnps))
    tmp = c(1, lABF)  # add on extra for null model
    my.denom = coloc:::logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    tmp1/sum(tmp1)
}

#' Posterior probabilities of causality from marginal Z-scores
#'
#' This function converts Z-scores to posterior probabilities of causality
#' i.e. not including the null model of no genetic effects,
#' so that the sum of the posterior probabilities for all variants is 1
#' @title Find posterior probabilities of causality
#' @param z Vector of marginal Z-scores
#' @param V Variance of the estimated effect size (can be obtained using Var.beta.cc function)
#' @param W Prior for the standard deviation of the effect size parameter, beta (W = 0.2 default)
#' @export
#' @return Vector of posterior probabilities
ppfunc <- function(z, V, W = 0.2) {
    r = W^2/(W^2 + V)
    bf = 0.5 * (log(1 - r) + (r * z^2))
    denom = coloc:::logsum(bf)
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp/sum(pp.tmp)
}

#' Posterior probabilities of causality from matrix of marginal Z-scores (1 simulation per row)
#'
#' This function converts a matrix of Z-scores (one row per simulation) to posterior probabilities of causality,
#' not including the null model of no genetic effects,
#' so that the sum of the posterior probabilities for each simulation (each row) is 1.
#' @title Find posterior probabilities of causality from matrix of Z-scores
#' @param zstar Matrix of marginal z-scores, one replicate per row
#' @param V Variance of the estimated effect size, one element per column of zstar
#' @param W Prior for the standard deviation of the effect size parameter, beta
#' @return Matrix of posterior probabilities of causality
#' @export
#' @author Chris Wallace
ppfunc.mat <- function(zstar, V, W = 0.2) {
    r = matrix(W^2/(W^2 + V), nrow = nrow(zstar), ncol = ncol(zstar), byrow = TRUE)  # see wakefield paper
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = apply(bf, 1, coloc:::logsum)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    pp.tmp = exp(bf - matrix(denom, nrow = nrow(bf), ncol = ncol(bf)))  # convert back from log scale
    pp.tmp/rowSums(pp.tmp)
}
