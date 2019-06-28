#' Simulate marginal z-scores (\eqn{Z_m}) from the joint z-scores (\eqn{Z_j}) using \eqn{E(Z_m) = Z_j \times \Sigma} and
#' \eqn{Z* \sim MVN(E(Z_m), \Sigma)}
#'
#' @title Simulate marginal from joint Z-score vector
#' @param Zj Vector of joint Z-scores (a vector of 0s except at the causal variant)
#' @param Sigma SNP correlation matrix
#' @param nrep Number of Z-score systems to simulate
#' @export
#' @return Matrix of simulated z-scores, one simulation per row
z_sim <- function(Zj, Sigma, nrep) {
    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    mexp.zm + ERR  # nrep is rows, nsnps is cols
}

#' Simulate nrep marginal Z-scores from joint Z-scores and convert these to posterior probabilities of causality
#'
#' Does not include posterior probabilities for null model so, the output will sum to 1
#' @title Simulate marginal from joint Z-score vector
#' @param Zj Vector of joint z-scores (0s except at CV)
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function)
#' @param nrep Number of marginal z-scores to simulate
#' @param W Prior for the standard deviation of the effect size parameter, beta
#' @param Sigma SNP correlation matrix
#' @export
#' @return Matrix of simulated posterior probabilties, one simulation per row
zj_pp <- function(Zj, V, nrep = 5000, W = 0.2, Sigma) {
    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    zstar = mexp.zm + ERR
    r = W^2/(W^2 + V)
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp/rowSums(pp.tmp)
}
