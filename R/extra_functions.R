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
    mean(x$covered)
}

#' Estimate for true effect at CV
#'
#' This function approximates the true effect at the causal variant given the expected value of the absolute marginal z-scores. For the correction we input the absolute marginal z-scores normalised by the corresponding posterior probabilities, \eqn{sum(abs(z0)*pp0)}.
#' @rdname mu.est
#' @title mu_est
#' @param X Some estimate of \eqn{E(|Z|)}
#' @export
#' @return Estimate of true effect at CV
#' @example
#' y = seq(0, 20, 0.005)
#' x <- sapply(y, function(m) mean(abs(rnorm(50000, mean = m))))
#' par(mfrow=c(1,2))
#' plot(x, y, xlab='E(|z|)', ylab='mu', main = 'Find mu value given expected |Z|')
#' abline(0, 1, col = 2)
#' plot(x, y, xlim=c(0,2), ylim=c(0,2), xlab='E(|z|)', ylab='mu', main = 'Zoom in to bottom left')
#' abline(0, 1, col = 2)
#' @author Anna Hutchinson
mu_est <- function(X) {
    y = seq(0, 20, 0.005)
    x <- sapply(y, function(m) mean(abs(stats::rnorm(50000, mean = m))))
    stats::approx(x, y, xout = X)$y
}
