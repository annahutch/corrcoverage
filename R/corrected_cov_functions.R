#' Obtain a credible set using the Bayesian approach for fine-mapping ([Maller et al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/23104008)).
#'
#' If the CV parameter is supplied, the output includes a binary indicator of whether the CV is contained in the set
#' @title credset
#' @param pp Vector of posterior probabilities
#' @param CV Optional parameter: Index of CV
#' @param thr Minimum threshold for credible set size
#' @export
#' @return  data.frame of claimed.cov (cumulative sum of the posterior probabilities of the variants forming the credible set), binary covered indicator (1 if CV is contained in the credible set), nvar (number of variants in the set)
credset <- function(pp, CV, thr) {
    o <- order(pp, decreasing = TRUE)  # order index for true pp
    cumpp <- cumsum(pp[o])  # cum sums of ordered pps
    wh <- which(cumpp > thr)[1]  # how many needed to exceed thr
    size <- cumpp[wh]
    if (missing(CV)) {
        data.frame(claimed.cov = size, nvar = wh)
    } else {
        contained = as.numeric(CV %in% o[1:wh])
        data.frame(claimed.cov = size, covered = contained, nvar = wh)
    }
}

#' Quicker credset function for matrix of posterior probabilities (using Cpp)
#'
#' @param pp Matrix of posterior probabilities (row per system)
#' @param CV Vector of indices of the causal variant for each row of posterior probabilities
#' @param thr Threshold
#'
#' @return Data.frame of claimed coverage (sum of posterior probabilities of variants in the set), binary covered indicator and number of variants.
#' @useDynLib corrcoverage
#' @importFrom Rcpp sourceCpp
#' @export
credset3 <- function(pp, CV=iCV, thr=0.6) {
  ret  <-  credsetmat(pp,CV,thr) ## list 1 = wh, 2 = size, 3=contained
  data.frame(claimed.cov=ret[[2]], covered=ret[[3]], nvar=ret[[1]])
}

#' Provide a corrected coverage estimate of the causal variant in the credible set
#'
#' Requires an estimate of the true effect at the CV (e.g. use maximum absolute z-score or output from mu_est function)
#' @rdname corrected_cov
#' @title corrected_cov
#' @param mu The true effect at the CV
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function)
#' @param Sigma SNP correlation matrix
#' @param pp0 Posterior probabilities of system of interest
#' @param thresh Minimum threshold for fine-mapping experiment
#' @param nrep Number of posterior probability systems to simulate for each variant considered causal
#' @export
#' @return Corrected coverage estimate
corrected_cov <- function(mu, V, Sigma, pp0, thresh, nrep = 1000) {
  nsnps <- length(pp0)

  # form joint z-score vectors
  temp <- diag(x = mu, nrow = nsnps, ncol = nsnps)
  zj <- lapply(seq_len(nrow(temp)), function(i) temp[i,]) # nsnp zj vectors for each snp considered causal

  # simulate pp systems
  pps <- mapply(zj_pp, zj, V, MoreArgs = list(nrep = nrep, Sigma = LD), SIMPLIFY = FALSE)

  # consider different CV as causal in each list
  n_pps <- length(pps)
  args <- 1:nsnps

  # obtain credible set for each simulation
  prop_cov <- lapply(1:n_pps, function(x) {
    tmp = credset_cov(pps[[x]], CV = rep(args[x], nrep), thr = thresh)
    sum(tmp)/length(tmp)
  }) %>% unlist()

  # final corrected coverage value
  sum(prop_cov * pp0)
}


#' Obtain corrected coverage estimate
#'
#' This function only requires the marginal summary statistics from GWAS
#' @rdname corrcov
#' @title corrcov
#' @param z Marginal z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param rep The number of simulated posterior probability systems to consider for the corrected coverage estimate. Simulations show that rep=1000 is a robust choice.
#' @export
#' @return Corrected coverage estimate
corrcov <- function(z, f, N0, N1, Sigma, thr, W = 0.2, rep = 1000) {
  ph0.tmp <- z0_pp(z, f, type = "cc", N = N0 + N1, s = N1/(N0 + N1), W = W)
  ph0 <- ph0.tmp[1]  # prob of the null
  pp0dash <- ph0.tmp[-1]  # pps of variants

  varbeta <- coloc:::Var.data.cc(f, N = N0 + N1, s = N1/(N0 + N1))  # variance of estimated effect size

  pp0 <- ppfunc(z, V = varbeta, W = W)  # posterior probs of system

  muhat <- est_mu(z, f, N0, N1)

  corrected_cov(mu = muhat, V = varbeta, Sigma = Sigma, pp0 = pp0, thresh = thr, nrep = rep)
}
