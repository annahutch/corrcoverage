#' Wakefield's lABF with prior SD as a parameter
#'
#' This function converts p-values to lABF
#' @title approx.bf.p
#' @param p p-values
#' @param f Minor allele frequencies
#' @param type Type of experiment ("quant" or "cc")
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=="quant"
#' @param W Prior for the standard deviation of the effect size parameter beta
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
#' @param W Prior for the standard deviation of the effect size parameter beta
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
#' @param W Prior for the standard deviation of the effect size parameter beta
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
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @export
#' @return Vector of posterior probabilities
ppfunc <- function(z, V, W=0.2) {
  r <- W^2/(W^2 + V)
  bf = 0.5 * (log(1 - r) + (r * z^2))
  denom <- coloc:::logsum(bf)
  pp.tmp <- exp(bf - denom) # convert back from log scale
  pp.tmp / sum(pp.tmp)
}

##' Bayesian finemapping analysis
##'
##' This function calculates posterior probabilities of different
##' causal variant for a single trait.
##'
##' If regression coefficients and variances are available, it
##' calculates Bayes factors for association at each SNP.  If only p
##' values are available, it uses an approximation that depends on the
##' SNP's MAF and ignores any uncertainty in imputation.  Regression
##' coefficients should be used if available.
##'
##' @title Bayesian finemapping analysis
##' @param dataset a list with the following elements
##' \describe{
##'
##'   \item{pvalues}{P-values for each SNP in dataset 1}
##'
##'   \item{N}{Number of samples in dataset 1}
##'
##'   \item{MAF}{minor allele frequency of the variants}
##'
##' \item{beta}{regression coefficient for each SNP from dataset 1}
##'
##' \item{varbeta}{variance of beta}
##'
##' \item{type}{the type of data in dataset 1 - either "quant" or "cc" to denote quantitative or case-control}
##'
##' \item{s}{for a case control dataset, the proportion of samples in dataset 1 that are cases}
##'
##'  \item{sdY}{for a quantitative trait, the population standard deviation of the trait.  if not given, it can be estimated from the vectors of varbeta and MAF}
##'
##' \item{snp}{a character vector of snp ids, optional. If present, it will be used to merge dataset1 and dataset2.  Otherwise, the function assumes dataset1 and dataset2 contain results for the same SNPs in the same order.}
##'
##' }
##'
##' Some of these items may be missing, but you must give
##' \itemize{
##' \item{always}{\code{type}}
##' \item{if \code{type}=="cc"}{\code{s}}
##' \item{if \code{type}=="quant" and \code{sdY} known}{\code{sdY}}
##' \item{if \code{type}=="quant" and \code{sdY} unknown}{\code{beta}, \code{varbeta}, \code{N}, \code{MAF}}
##' and then either
##' \item{}{\code{pvalues}, \code{MAF}}
##' \item{}{\code{beta}, \code{varbeta}}
##' }
##'
##'
##' @param p1 prior probability a SNP is associated with the trait 1, default 1e-4
##' @return a \code{data.frame}:
##' \itemize{
##' \item an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability of the SNP being causal
##' }
##' @author Chris Wallace
##' @export
finemap.abf <- function(dataset, p1=1e-4) {

  if(!is.list(dataset))
    stop("dataset must be a list.")

  df <- process.dataset(d=dataset, suffix="")
  nsnps <- nrow(df)
  dfnull <- df[1,]
  for(nm in colnames(df))
    dfnull[,nm] <- NA
  dfnull[,"snp"] <- "null"
  dfnull[,"lABF."] <- 1
  df <- rbind(df,dfnull)
  ## data.frame("V."=NA,
  ##            z.=NA,
  ##            r.=NA,
  ##            lABF.=1,
  ##            snp="null"))
  df$prior <- c(rep(p1,nsnps),1-nsnps*p1)

  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- logsum(df$lABF + df$prior)
  df$SNP.PP <- exp(df$lABF - my.denom.log.abf)

  return(df)
}

#' Obtain posterior probabilities for each SNP causal from multiple marginal z-scores simulations
#'
#' This function converts a matrix of z-scores (one row per simulation) to posterior probabilities, not including the null model of no genetic effects, so that the sum of the posterior probabilities for each simulation (each row) is 1.
#' @title ppfunc.mat
#' @param zstar Matrix of marginal z-scores, one replicate per row
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function), one element per column of zstar
#' @param W Prior for the standard deviation of the effect size parameter beta
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
