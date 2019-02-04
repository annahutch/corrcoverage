#' Variance of the estimated effect size, V=var(\hat{\beta}), for case-control
#'
#' @title Var.data.cc
#' @param f Minor allele frequencies
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1)
#' @return Variance of estimated effect size, V.
#' @author Claudia Giambartolomei
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
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
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

#' Correlation matrix of SNPS
#'
#' @title cor2
#' @param x Phased haplotype matrix, rows as samples and columns as SNPs
#' @return Correlation matrix
#' @author Chris Wallace
cor2 <- function (x) {
  1/(NROW(x) - 1) * crossprod(scale(x, TRUE, TRUE))
}

#' Obtain phased haplotype matrix from UK10K data
#'
#' @title geth
#' @param nsnp Number of SNPs to consider for fine-mapping
#' @return Haplotype matrix
#' @author Chris Wallace
geth <- function(nsnp) {
  ## real data from UK10K
  file.ldd="/home/cew54/share/Data/reference/lddetect/EUR/fourier_ls-chr22.bed"
  file.vcf="/home/cew54/share/Data/reference/UK10K/BCF/chr22.bcf.gz"

  ## ldblocks
  ldd <- fread(file.ldd)

  ## split bcf by ldblocks
  ldd[,blocknum:=1:.N]
  ldd[,dist:=stop-start]
  ldd[,comm:=paste0("/home/cew54/share/bin/bcftools view ",file.vcf,
                    " --min-af 0.02:minor --max-alleles 2 --min-alleles 2 ",
                    " -r chr",22,":",start,"-",stop," -Ov ")] # -o ",tmp)]
  gethap <- function(i) {
    y=fread(ldd$comm[i])
    ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
    rownames(ha) <- paste0("pos",y$POS)
    t(ha)
  }

  block <- sample(which(ldd$dist < 1.2e+6),1) # use smallest LD block to be fast
  h <- gethap(block) # rows=samples, cols=snps
  use <- apply(h,2,var)>0 & colMeans(h) > 0.01 & colMeans(h)<0.99  # no monomorphs
  h <- h[,use,drop=FALSE]
  nmax <-floor(nsnp/2) # keep this small to make simulations fast
  if(ncol(h)>nmax) {
    start <- sample(1:ncol(h), 1) # choose a random starting point
    start2 <- start+floor(ncol(h)/2) # choose another random starting point, the other end
    idx<-c(start+1:nmax, start2+1:nmax) # doing this means we're not just looking at one LD block of SNPs, so there will be some snps completely seperate to CV
    idx<-idx %% ncol(h)
    h <- h[,idx]
  }
  list("h"=h,"block"=block)
}

#' Wakefield's ABF with prior SD as a parameter
#'
#' @title approx.bf.p
#' @param p p-values
#' @param f Minor allele frequencies
#' @param type Type of experiment ("quant" or "cc")
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=="quant"
#' @param W Prior for the standard deviation of the effect size parameter $\beta$
#' @return data.frame containing lABF and intermediate calculations
#' @author Chris Wallace
approx.bf.p <- function (p, f, type, N, s, W) {
  sd.prior <- W
  V <- coloc::Var.data.cc(f, N, s)
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  return(ret)
}

#' Quick posterior probabilities from p-values with prior SD as a parameter
#'
#' This derivation also includes the null model of no genetic effects
#' @title p.vals_pp
#' @param pvals p-values of SNPs
#' @param f Minor allele frequencies
#' @param type Type of experiment ("quant" or "cc")
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=="quant"
#' @param W Prior for the standard deviation of the effect size parameter $\beta$
#' @return Posterior probabilities of null model (no genetic effect) and causality for each SNP
#' @author Anna Hutchinson
p.vals_pp <- function(pvals, f, type, N, s, W=0.2){
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
#' This derivation also includes the null model of no genetic effects
#' @title z0_pp
#' @param z0 Marginal z-scores of SNPs
#' @param f Minor allele frequencies
#' @param type Type of experiment ("quant" or "cc")
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1), ignored if type=="quant"
#' @param W Prior for the standard deviation of the effect size parameter $\beta$
#' @return Posterior probabilities of null model (no genetic effect) and causality for each SNP
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
#' This derivation does not include the null model, so one posterior probability is obtained for each marginal
#' z-score and these sum to 1.
#' @title Posterior probabilities for vectors
#' @param z Vector of marginal z-scores
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function)
#' @param W Prior for the standard deviation of the effect size parameter $\beta$
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
#' Each row should be a simulation and each column a SNP
#'
#' @title Posterior probabilities for matrices
#' @param zstar Matrix of marginal z-scores, one replicate per row
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function), one element per column of zstar
#' @param W Prior for the standard deviation of the effect size parameter $\beta$
#' @return Vector of posterior probabilities
#' @author Chris Wallace
ppfunc.mat <- function(zstar, V, W=0.2){
  r <- matrix(omega^2/(omega^2 + V), nrow=nrow(zstar), ncol=ncol(zstar), byrow=TRUE) # see wakefield paper
  bf = 0.5 * (log(1 - r) + (r * zstar^2))
  denom <- apply(bf, 1, coloc:::logsum) # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
  pp.tmp <- exp(bf - matrix(denom, nrow=nrow(bf),ncol=ncol(bf))) # convert back from log scale
  pp.tmp / rowSums(pp.tmp)
}
