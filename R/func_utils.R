#' Internal function: Simulate nrep ABFs from joint Z-score vector
#'
#' Does not include posterior probabilities for null model
#' @param Zj joint z vector
#' @param int.Sigma internal sigma
#' @param int.nrep internal nrep
#' @param int.ERR internal error matrix
#' @param int.r internal r
#' @return Matrix of simulated ABFs, one simulation per row
.zj_abf = function(Zj, int.Sigma, int.nrep, int.ERR, int.r) {
  
  stopifnot(inherits(Zj, "numeric"))
  stopifnot(inherits(int.ERR, "matrix"))
  stopifnot(inherits(int.r, "numeric"))
  exp.zm = Zj %*% int.Sigma
  mexp.zm = matrix(exp.zm, int.nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
  zstar = mexp.zm + int.ERR
  0.5 * t(log(1 - int.r) + (int.r * t(zstar^2)))
}

#' Internal function: Simulate nrep posterior probabilities of causality from joint Z-score vector
#'
#' Does not include posterior probabilities for null model
#' @title Simulate posterior probabilities of causality from joint
#'     Z-score vector
#' @param Zj joint z vector
#' @param int.Sigma internal sigma
#' @param int.nrep internal nrep
#' @param int.ERR internal error matrix
#' @param int.r internal r
#' @return Matrix of simulated posterior probabilties of causality,
#'     one simulation per row
.zj_pp = function(Zj, int.Sigma, int.nrep, int.ERR, int.r) {
  stopifnot(inherits(Zj, "numeric"))
  stopifnot(inherits(int.ERR, "matrix"))
  stopifnot(inherits(int.r, "numeric"))
  zj_pp_arma(Zj, int.Sigma, int.nrep, int.ERR, int.r)
  # exp.zm = Zj %*% int.Sigma
  # mexp.zm = matrix(exp.zm, int.nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
  # zstar = mexp.zm + int.ERR
  # bf = 0.5 * t(log(1 - int.r) + (int.r * t(zstar^2)))
  # denom = logsum_matrix(bf) # faster
  # exp(bf - denom)  # convert back from log scale
}

#' @importFrom matrixStats rowMaxs
NULL

##' matrix-ified version of logsum to avoid needing apply()
##'
##' @title logsum rows of a matrix
##' @param x numeric matrix
##' @return rowwise sums
##' @author Chris Wallace
logsum_matrix=function(x) {
  my.max=rowMaxs(x)
  my.max + log(rowSums(exp(x - my.max)))
}

