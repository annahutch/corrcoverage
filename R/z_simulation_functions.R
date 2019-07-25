#' Simulate marginal z-scores (\eqn{Z_m}) from the joint z-scores (\eqn{Z_j}) using \eqn{E(Z_m) = Z_j \times \Sigma} and
#' \eqn{Z* \sim MVN(E(Z_m), \Sigma)}
#'
#' @title Simulate marginal Z-scores from joint Z-score vector
#' @param Zj Vector of joint Z-scores (a vector of 0s except at the CV)
#' @param Sigma SNP correlation matrix
#' @param nrep Number of Z-score systems to simulate
#' @return Matrix of simulated posterior probabilties, one simulation per row
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#'
#' # derive joint Z score vector
#' Zj <- rep(0, nsnps)
#' iCV <- 4 # index of CV
#' mu <- 5 # true effect at CV
#' Zj[iCV] <- mu
#'
#' ## generate example LD matrix
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
#' LD <- cor2(X)
#'
#' res <- z_sim(Zj, Sigma = LD, nrep = 100)
#' res[c(1:5), c(1:5)]
#'
#' @export
#' @author Anna Hutchinson
z_sim <- function(Zj, Sigma, nrep) {
    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    mexp.zm + ERR  # nrep is rows, nsnps is cols
}

#' Simulate nrep marginal Z-scores from joint Z-scores and convert these to posterior probabilities of causality
#'
#' Does not include posterior probabilities for null model
#' @title Simulate posterior probabilities of causality from joint Z-score vector
#' @param Zj Vector of joint Z-scores (0s except at CV)
#' @param V Variance of the estimated effect size (can be obtained using Var.beta.cc function)
#' @param nrep Number of posterior probability systems to simulate (default 1000)
#' @param W Prior for the standard deviation of the effect size parameter, beta (default 0.2)
#' @param Sigma SNP correlation matrix
#' @return Matrix of simulated posterior probabilties, one simulation per row
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#' Zj <- rep(0, nsnps)
#' iCV <- 4 # index of CV
#' mu <- 5 # true effect at CV
#' Zj[iCV] <- mu
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
#' LD <- cor2(X)
#' maf <- colMeans(X)
#'
#' ## generate V (variance of estimated effect sizes)
#' varbeta <- Var.data.cc(f = maf, N = 5000, s = 0.5)
#'
#' res <- zj_pp(Zj, V = varbeta, nrep = 5, W = 0.2, Sigma = LD)
#'
#' res[c(1:5), c(1:5)]
#'
#'
#' @export
#' @author Anna Hutchinson
zj_pp <- function(Zj, V, nrep = 1000, W = 0.2, Sigma) {
    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    zstar = mexp.zm + ERR
    r = W^2/(W^2 + V)
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = apply(bf,1,logsum)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    exp(bf - denom)  # convert back from log scale
}

#' Internal function: Simulate nrep posterior probabilities of causality from joint Z-score vector
#'
#' Does not include posterior probabilities for null model
#' @title Simulate posterior probabilities of causality from joint Z-score vector
#' @param Zj Vector of joint Z-scores (0s except at CV)
#' @return Matrix of simulated posterior probabilties of causality, one simulation per row
pp_ERR = function(Zj) {
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    zstar = mexp.zm + ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = apply(bf, 1, logsum) # get different denom for each rep
    exp(bf - denom)  # convert back from log scale
}
