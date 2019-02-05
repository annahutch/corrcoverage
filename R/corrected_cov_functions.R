#' Obtain a credible set using the Bayesian approach for fine-mapping ([Maller et al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/23104008)).
#'
#' If the CV parameter is supplied, the output includes a binary indicator of whether the CV is contained in the set
#' @title credset
#' @param pp Vector of posterior probabilities
#' @param CV Index of CV
#' @param thr Minimum threshold for credible set size
#' @return  data.frame of claimed.cov (cumulative sum of the posterior probabilities of the variants forming the credible set), binary covered indicator (1 if CV is contained in the credible set), nvar (number of variants in the set)
credset <- function(pp, CV = iCV, thr) {
  o <- order(pp, decreasing = TRUE)  # order index for true pp
  cumpp <- cumsum(pp[o])  # cum sums of ordered pps
  wh <- which(cumpp > thr)[1]  # how many needed to exceed thr
  size <- cumpp[wh]
  if(missing(CV)) {
    data.frame(claimed.cov = size, nvar = wh)
  } else{
    contained=as.numeric(CV %in% o[1:wh])
    data.frame(claimed.cov = size, covered = contained, nvar = wh)
  }
}


#' Provide a corrected coverage estimate of the causal variant in the credible set
#'
#' Requires an estimate of the true effect at the CV
#' @rdname corrected_cov
#' @title corrected_cov
#' @param mu The true effect at the CV
#' @param nsnps Number of SNPs
#' @param V Variance of the estimated effect size (can be obtained using var.beta.cc function)
#' @param Sigma SNP correlation matrix
#' @param pp0 Posterior probabilities of system of interest
#' @param thresh Minimum threshold for fine-mapping experiment
#' @return Corrected coverage estimate
corrected_cov <- function(mu, nsnps = 200, V, Sigma, pp0, thresh) {
  temp <- diag(x = mu, nrow = nsnps, ncol = nsnps)
  zj <- do.call(c, apply(temp, 1, list))  # nsnp zj vectors for each snp considered causal

  # simulate pp systems
  pps <- mapply(zj_pp, zj, V, MoreArgs = list(Sigma = LD), SIMPLIFY = FALSE)

  # consider different CV as causal in each list
  n_pps <- length(pps)
  args <- 1:nsnps

  # obtain credible set for each simulation
  d5 <- lapply(1:n_pps, function(x) {
    apply(pps[[x]], 1, credset, args[x], thr = thresh) %>% data.table:::rbindlist()
  })

  invlogit <- function(x) exp(x)/(1 + exp(x))
  logit <- function(x) log(x/(1 - x))

  # resize claimed coverage resize so don't have 0 on denominator
  resize <- function(x) {
    x[x > 0.99999999] <- 0.999999
    return(x)
  }

  claim.cov <- lapply(d5, function(p) resize(p$claimed.cov))
  logitclaimed <- lapply(claim.cov, function(p) logit(p))
  y <- mapply(cbind, d5, logit.claim = logitclaimed, SIMPLIFY = FALSE)

  # if get error when fitting model use mean(covered)
  model <- function(y) {
    out <- tryCatch({
      lapply(y, pred_logit) %>% unlist()
    }, error = function(cond) {
      lapply(y, pred_na) %>% unlist()
    })
    return(out)
  }

  final <- model(y)

  # if NaNs, fit intercept only model

  if (mean(final) == "NaN") {
    final1 <- lapply(d5, pred_na) %>% unlist()
  } else {
    final1 <- final
  }

  # final corrected coverage value
  sum(final1 * pp0)
}


#' Obtain corrected coverage estimate
#'
#' This function only requires the marginal summary statistics from GWAS
#' @rdname corrcov
#' @title corrcov
#' @param z0 Marginal z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment
#' @return Corrected coverage estimate
corrcov <- function(z0, f, N0, N1, Sigma, thr) {
  ph0.tmp <- z0_pp(z0, f, type = "cc", N = N0+N1, s = N1/(N0+N1), W = 0.2)
  ph0 <- ph0.tmp[1]  # prob of the null
  pp0dash <- ph0.tmp[-1]  # pps including the null

  varbeta <- Var.data.cc(f, N0+N1, N1/(N0+N1))  # variance of beta

  pp0 <- ppfunc(z0, V = varbeta)  # posterior probs of system

  muhat.gam <- mu_est(sum(abs(z0) * pp0))  # estimate for true effect at CV

}
