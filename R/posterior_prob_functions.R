#' Wakefield's ABF with prior SD as a parameter
#'
#' This function converts p-values to ABF
#' @title approx.bf.p
#' @param p p-values
#' @param f Minor allele frequencies
#' @param type Type of experiment ("quant" or "cc")
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=="quant"
#' @param W Prior for the standard deviation of the effect size parameter \beta
#' @return data.frame containing lABF and intermediate calculations
#' @export
#' @author Chris Wallace
approx.bf.p <- function (p, f, type, N, s, W) {
  sd.prior <- W
  V <- coloc:::Var.data.cc(f, N, s)
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  return(ret)
}

#' Quick posterior probabilities from p-values with prior SD as a parameter
#'
#' This function converts p-values to posterior probabilities, including the null model of no genetic effects
#' @title pvals_pp
#' @param pvals p-values of SNPs
#' @param f Minor allele frequencies
#' @param type Type of experiment ("quant" or "cc")
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=="quant"
#' @param W Prior for the standard deviation of the effect size parameter \beta
#' @return Posterior probabilities of null model (no genetic effect) and causality for each SNP
#' @export
#' @author Anna Hutchinson
pvals_pp <- function(pvals, f, type, N, s, W=0.2){
  tmp <- approx.bf.p(p = pvals, f = f, type = "cc", N = N, s = 0.5, W = W)[,"lABF"]
  p1 <- 1e-04 # hard coded - bad code
  nsnps <- length(tmp)
  prior <- c(1-nsnps*p1,rep(p1,nsnps))
  tmp <- c(1,tmp) # add on extra for null model
  my.denom <- coloc:::logsum(tmp + prior)
  tmp1 <- exp(tmp+prior - my.denom)
  tmp1/sum(tmp1)
}

#' Quick posterior probabilities from marginal z-scores with prior SD as a parameter
#'
#' This function converts z-scores to posterior probabilities, including the null model of no genetic effects
#' @title z0_pp
#' @param z0 Marginal z-scores of SNPs
#' @param f Minor allele frequencies
#' @param type Type of experiment ("quant" or "cc")
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=="quant"
#' @param W Prior for the standard deviation of the effect size parameter \beta
#' @return Posterior probabilities of null model (no genetic effect) and causality for each SNP
#' @export
#' @author Anna Hutchinson
z0_pp <- function(z0, f, type, N, s, W = 0.2){
  pvals <- pnorm(abs(z0),lower.tail = FALSE)*2 # convert z-scores to p-values
  tmp <- approx.bf.p(p = pvals, # find approx bf
                     f = f, type = "cc", N = N, s = 0.5, W = W)[,"lABF"]
  p1 <- 1e-04 # hard coded - bad code
  nsnps <- length(tmp)
  prior <- c(1-nsnps*p1,rep(p1,nsnps))
  tmp <- c(1,tmp) # add on extra for null model
  my.denom <- coloc:::logsum(tmp + prior)
  tmp1 <- exp(tmp+prior - my.denom)
  tmp1/sum(tmp1)
}

#' Quick posterior probabilities for each SNP causal from marginal z-scores
#'
#' This function converts z-scores to posterior probabilities, not including the null model of no genetic effects, so that the sum of the posterior probabilities is 1
#' @title ppfunc
#' @param z Vector of marginal z-scores
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function)
#' @param W Prior for the standard deviation of the effect size parameter \beta
#' @export
#' @return Vector of posterior probabilities
ppfunc <- function(z, V, W=0.2) {
  r <- W^2/(W^2 + V)
  bf = 0.5 * (log(1 - r) + (r * z^2))
  denom <- coloc:::logsum(bf)
  pp.tmp <- exp(bf - denom) # convert back from log scale
  pp.tmp / sum(pp.tmp)
}

#' Obtain posterior probabilities for each SNP causal from multiple marginal z-scores simulations
#'
#' This function converts a matrix of z-scores (one row per simulation) to posterior probabilities, not including the null model of no genetic effects, so that the sum of the posterior probabilities for each simulation (each row) is 1.
#' @title ppfunc.mat
#' @param zstar Matrix of marginal z-scores, one replicate per row
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function), one element per column of zstar
#' @param W Prior for the standard deviation of the effect size parameter \beta
#' @return Vector of posterior probabilities
#' @export
#' @author Chris Wallace
ppfunc.mat <- function(zstar, V, W=0.2){
  r <- matrix(W^2/(W^2 + V), nrow=nrow(zstar), ncol=ncol(zstar), byrow=TRUE) # see wakefield paper
  bf = 0.5 * (log(1 - r) + (r * zstar^2))
  denom <- apply(bf, 1, coloc:::logsum) # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
  pp.tmp <- exp(bf - matrix(denom, nrow=nrow(bf),ncol=ncol(bf))) # convert back from log scale
  pp.tmp / rowSums(pp.tmp)
}
