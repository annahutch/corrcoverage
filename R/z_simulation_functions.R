#' Simulate marginal z-scores from the joint z-scores using \eqn{Z_m ~ MVN(E(Z_j),\Sigma)}
#'
#' @title z_sim
#' @param Zj Vector of joint z-scores (0s except at CV)
#' @param Sigma SNP correlation matrix
#' @param nrep Number of simulated z-scores
#' @return Matrix of simulated z-scores, one simulation per row
z_sim <- function(Zj, Sigma, nrep) {
    exp.zm = Zj %*% Sigma  # find E(X_m) (of for each SNP being causal)
    mvtnorm:::rmvnorm(nrep, exp.zm, Sigma)  # nrep is rows, nsnps is cols
}

#' Simulate nrep marginal z-scores from joint z-scores and convert these to posterior probabilities
#'
#' Does not include posterior probabilities for null model
#' @title zj_pp
#' @param Zj Vector of joint z-scores (0s except at CV)
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function)
#' @param nrep Number of marginal z-scores to simulate
#' @param W Prior for the standard deviation of the effect size parameter $\beta$
#' @param Sigma SNP correlation matrix
#' @return Matrix of simulated posterior probabilties, one simulation per row
zj_pp <- function(Zj, V, nrep = 5000, W = 0.2, Sigma) {
    exp.zm = Zj %*% Sigma  # find E(Z_m)
    zstar = mvtnorm:::rmvnorm(nrep, exp.zm, Sigma)  # nrep is rows, nsnps is cols
    r <- W^2/(W^2 + V)
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom <- coloc:::logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    pp.tmp <- exp(bf - denom)  # convert back from log scale
    pp.tmp/rowSums(pp.tmp)
}
