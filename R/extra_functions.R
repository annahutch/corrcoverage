#' Variance of the estimated effect size for case-control data
#'
#' @title Var.data.cc
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
#' @title cor2
#' @param x Phased haplotype matrix, rows as samples and columns as SNPs
#' @return Correlation matrix
#' @export
#' @author Chris Wallace
cor2 <- function(x) {
    1/(NROW(x) - 1) * crossprod(scale(x, TRUE, TRUE))
}

#' Build and predict corrected coverage using logistic GAM
#'
#' This function builds a GAM model with log(claimed_coverage/(1-claimed_coverage)) (logit(claimed) as the predictor and the binary covered value as the response. This model is then used to accurately predict the coverage probability of the causal variant in the credible set
#' @rdname pred_logit
#' @title pred_logit
#' @param x data.frame with column for 'covered' and 'logit.claim'
#' @param size size of the credible set to correct (sum of the posterior probabilities of the variants)
#' @export
#' @return Predicted probability of covered
#' @author Anna Hutchinson
pred_logit <- function(x, size) {
    m <- mgcv::gam(covered ~ s(logit.claim), data = x, family = "binomial")
    invlogit(predict(m, newdata = data.frame(logit.claim = logit(size))))
}

#' Prediction for corrected coverage if gam cannot be fitted
#'
#' Predicted corrected coverage is the proportion of simulated credible sets that contain the causal variant
#' @rdname pred_na
#' @title pred_na
#' @param x data.frame with a binary 'covered' column
#' @return Predicted coverage
#' @export
#' @author Anna Hutchinson
pred_na <- function(x) {
  sum(x$covered)/length(x$covered)
}

#' Find an estimate for the true effect at the causal variant
#'
#' @param z Vector of Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of cases
#' @param N1 Number of controls
#'
#' @return Estimate of the true effect at the causal variant
#' @export
#'
#' @author Anna Hutchinson
est_mu <- function(z, f, N0, N1){
  ph0.tmp <- z0_pp(z = z, f = f, type = "cc", N = N0+N1, s = 0.5)
  ph0 <- ph0.tmp[1] # prob of the null
  mean(c(sum(abs(z)*ph0.tmp[-1]),(1-ph0.tmp[1])*max(abs(z))))
}
