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

#' Build and predict corrected coverage using logistic GAM
#'
#' This function builds a GAM model with log(claimed_coverage/(1-claimed_coverage)) (logit(claimed))
#' as the predictor and the binary covered value as the response.
#' This model is then used to predict the coverage probability of the causal variant in the credible set
#' @rdname pred_logit
#' @title pred_logit
#' @param x data.frame with column for 'covered' and 'logit.claim'
#' @param size size of the credible set to correct (sum of the posterior probabilities of causality of the variants)
#' @return Predicted probability of covered
#' @author Anna Hutchinson
pred_logit <- function(x, size) {
    m = mgcv::gam(covered ~ s(logit.claim), data = x, family = "binomial")
    invlogit(predict(m, newdata = data.frame(logit.claim = logit(size))))
}

#' Proportion of simulated credible sets containing the causal variant
#'
#' @rdname prop_cov
#' @title Proportion of credible sets containing the causal variant
#' @param x data.frame with a binary 'covered' column
#' @return Proportion of x with x = 1
#' @export
#' @author Anna Hutchinson
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
    p1 <- 1e-04  # hard coded
    nsnps <- length(lABF)
    prior <- c(1 - nsnps * p1, rep(p1, nsnps))
    tmp <- c(1, lABF)  # add on extra for null model
    my.denom <- coloc:::logsum(tmp + prior)
    tmp1 <- exp(tmp + prior - my.denom)
    ph0.tmp <- tmp1/sum(tmp1)

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
    my.denom = coloc:::logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    ph0.tmp = tmp1/sum(tmp1)
    ph0 = ph0.tmp[1]  # prob of the null
    mean(c(sum(abs(z) * ph0.tmp[-1]), (1 - ph0.tmp[1]) * max(abs(z))))
}
