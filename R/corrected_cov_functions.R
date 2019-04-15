#' Obtain a corrected coverage estimate of the causal variant in the credible set
#'
#' Requires an estimate of the true effect at the CV (e.g. use maximum absolute z-score or output from corrcoverage::mu_est function)
#' @rdname corrected_cov
#' @title Obtain a corrected coverage estimate of the causal variant in the credible set
#' @param mu The true effect at the CV
#' @param V Variance of the estimated effect size (can be obtained using Var.beta.cc function)
#' @param W Prior for the standard deviation of the effect size parameter beta (W=0.2 default)
#' @param Sigma SNP correlation matrix
#' @param pp0 Posterior probabilities of system of interest
#' @param thresh Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param nrep Number of posterior probability systems to simulate for each variant considered causal (nrep = 1000 default)
#' @export
#' @return Corrected coverage estimate
corrected_cov <- function(mu, V, W = 0.2, Sigma, pp0, thresh = 0.95, nrep = 1000) {

    nsnps = length(pp0)

    # form joint z-score vectors
    temp = diag(x = mu, nrow = nsnps, ncol = nsnps)
    zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal

    # simulate ERR matrix
    ERR = mvtnorm:::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    r = W^2/(W^2 + V)

    pp_ERR = function(Zj) {
        exp.zm = Zj %*% Sigma
        mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
        zstar = mexp.zm + ERR
        bf = 0.5 * (log(1 - r) + (r * zstar^2))
        denom = coloc:::logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
        pp.tmp = exp(bf - denom)  # convert back from log scale
        pp.tmp/rowSums(pp.tmp)
    }

    # simulate pp systems
    pps <- mapply(pp_ERR, zj, SIMPLIFY = FALSE)

    # consider different CV as causal in each list
    n_pps <- length(pps)
    args <- 1:nsnps

    # obtain credible set for each simulation
    d5 <- lapply(1:n_pps, function(x) {
        credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thresh)
    })

    prop_cov <- lapply(d5, prop_cov) %>% unlist()

    sum(prop_cov * pp0)
}


#' Obtain corrected coverage estimate using Z-scores and mafs
#'
#' This function only requires the marginal summary statistics from GWAS
#' @rdname corrcov
#' @title Obtain corrected coverage estimate using Z-scores and mafs
#' @param z Marginal Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 1000 default)
#' @export
#' @return Corrected coverage estimate
corrcov <- function(z, f, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000) {
    varbeta = 1/(2 * (N0 + N1) * f * (1 - f) * (N1/(N0 + N1)) * (1 - (N1/(N0 + N1))))
    r = W^2/(W^2 + varbeta)
    bf = 0.5 * (log(1 - r) + (r * z^2))
    p1 = 1e-04  # hard code
    nsnps = length(bf)
    prior = c(1 - nsnps * p1, rep(p1, nsnps))
    tmp = c(1, bf)  # add on extra for null model
    my.denom = coloc:::logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    ph0.tmp = tmp1/sum(tmp1)

    ph0 = ph0.tmp[1]  # prob of the null
    pp0dash = ph0.tmp[-1]  # pps of variants

    # posterior probs of true system
    pp.tmp = exp(bf - coloc:::logsum(bf))
    pp0 = pp.tmp/sum(pp.tmp)

    # estimate mu
    muhat = mean(c(sum(abs(z) * pp0dash), (1 - ph0) * max(abs(z))))

    #### corrected coverage
    temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
    zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal
    # simulate ERR matrix

    ERR = mvtnorm:::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    pp_ERR = function(Zj) {
        exp.zm = Zj %*% Sigma
        mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
        zstar = mexp.zm + ERR
        bf = 0.5 * (log(1 - r) + (r * zstar^2))
        denom = coloc:::logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
        pp.tmp = exp(bf - denom)  # convert back from log scale
        pp.tmp/rowSums(pp.tmp)
    }
    # simulate pp systems
    pps <- mapply(pp_ERR, zj, SIMPLIFY = FALSE)
    # consider different CV as causal in each list
    n_pps <- length(pps)
    args <- 1:nsnps

    # obtain credible set for each simulation
    d5 <- lapply(1:n_pps, function(x) {
        credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })

    prop_cov <- lapply(d5, prop_cov) %>% unlist()

    sum(prop_cov * pp0)
}

#' Obtain corrected coverage estimate using estimated effect sizes and their standard errors
#'
#' This function only requires the marginal summary statistics from GWAS
#' @rdname corrcov_bhat
#' @title Obtain corrected coverage estimate using estimated effect sizes and their standard errors
#' @param bhat Estimated effect sizes from single-SNP logistic regressions
#' @param V Variance of estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 1000 default)
#' @return Corrected coverage estimate
#' @export
#'
corrcov_bhat <- function(bhat, V, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000) {
    z <- bhat/sqrt(V)
    r <- W^2/(W^2 + V)
    bf = 0.5 * (log(1 - r) + (r * z^2))
    p1 = 1e-04  # hard code
    nsnps = length(bf)
    prior = c(1 - nsnps * p1, rep(p1, nsnps))
    tmp = c(1, bf)  # add on extra for null model
    my.denom = coloc:::logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    ph0.tmp = tmp1/sum(tmp1)

    ph0 = ph0.tmp[1]  # prob of the null
    pp0dash = ph0.tmp[-1]  # pps of variants

    # posterior probs of true system
    pp.tmp = exp(bf - coloc:::logsum(bf))
    pp0 = pp.tmp/sum(pp.tmp)

    # estimate mu
    muhat = mean(c(sum(abs(z) * pp0dash), (1 - ph0) * max(abs(z))))

    #### corrected coverage
    temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
    zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal
    # simulate ERR matrix

    ERR = mvtnorm:::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    pp_ERR = function(Zj) {
        exp.zm = Zj %*% Sigma
        mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
        zstar = mexp.zm + ERR
        bf = 0.5 * (log(1 - r) + (r * zstar^2))
        denom = coloc:::logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
        pp.tmp = exp(bf - denom)  # convert back from log scale
        pp.tmp/rowSums(pp.tmp)
    }

    # simulate pp systems
    pps <- mapply(pp_ERR, zj, SIMPLIFY = FALSE)
    # consider different CV as causal in each list
    n_pps <- length(pps)
    args <- 1:nsnps

    # obtain credible set for each simulation
    d5 <- lapply(1:n_pps, function(x) {
        credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })

    prop_cov <- lapply(d5, prop_cov) %>% unlist()

    sum(prop_cov * pp0)
}

#' @title Use simulated pps to find corrected coverage estimate
#'
#' @rdname quick_corrcov
#' @param thr Threshold value to exceed (default is 0.95)
#' @param simulated.pps A list of matrices of simulated posterior probabilities (1 simulation per row) where each list element is for the corresponding SNP considered causal
#' @param pp The posterior probabilities of the original system
#'
#' @return Corrected coverage estimate
#' @export
quick_corrcov <- function(thr = 0.95, simulated.pps, pp) {
    n_pps <- length(simulated.pps)
    args <- 1:nsnps

    d5 <- lapply(1:n_pps, function(x) {
        credsetC(simulated.pps[[x]], CV = rep(args[x], dim(simulated.pps[[x]])[1]), thr = thr)
    })

    prop_cov <- lapply(d5, prop_cov) %>% unlist()
    sum(prop_cov * pp)
}

#' @title Use simulated pps to find corrected coverage estimate and cred set
#'
#' @rdname quick_corrcov_cs
#' @param thr Threshold value to exceed (default is 0.95)
#' @param simulated.pps A list of matrices of simulated posterior probabilities (1 simulation per row) where each list element is for the corresponding SNP considered causal
#' @param pp The posterior probabilities of the original system
#'
#' @return A list of the credible set obtained using the specified threshold, the corrected coverage estimate of this credible set, the threshold the user specified and the size of the credible set (the sum of the pps of the variants)
#' @export
#'
quick_corrcov_cs <- function(thr = 0.95, simulated.pps, pp) {
    n_pps <- length(simulated.pps)
    args <- 1:nsnps

    d5 <- lapply(1:n_pps, function(x) {
        credsetC(simulated.pps[[x]], CV = rep(args[x], dim(simulated.pps[[x]])[1]), thr = thr)
    })

    prop_cov <- lapply(d5, prop_cov) %>% unlist()
    corr_cov <- sum(prop_cov * pp)

    o <- order(pp, decreasing = TRUE)
    cumpp <- cumsum(pp[o])
    wh <- which(cumpp > thr)[1]
    list(credset = names(pp)[o[1:wh]], corr.cov = corr_cov, thr = thr, size = cumpp[wh])
}

#' Obtain corrected coverage estimate using Z-scores and mafs (limiting simulations used for estimation to those with correct nvar)
#'
#' This function requires the marginal summary statistics from GWAS and an nvar value. It should only be used when nvar is very low ($<3$) and there is some evidence to suggest that only simulated credible sets with this nvar value should be used to derive the corrected coverage estimate.
#' @rdname corrcov_nvar
#' @title Obtain corrected coverage estimate using Z-scores and mafs
#' @param z Marginal Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param nvar The number of variants that simulated credible sets used for estimation should contain
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 10000 default because we trim out the ones without correct nvar)
#' @export
#' @return Corrected coverage estimate
corrcov_nvar <- function(z, f, N0, N1, Sigma, nvar, thr = 0.95, W = 0.2, nrep = 10000) {
  varbeta = 1/(2 * (N0 + N1) * f * (1 - f) * (N1/(N0 + N1)) * (1 - (N1/(N0 + N1))))
  r = W^2/(W^2 + varbeta)
  bf = 0.5 * (log(1 - r) + (r * z^2))
  p1 = 1e-04  # hard code
  nsnps = length(bf)
  prior = c(1 - nsnps * p1, rep(p1, nsnps))
  tmp = c(1, bf)  # add on extra for null model
  my.denom = coloc:::logsum(tmp + prior)
  tmp1 = exp(tmp + prior - my.denom)
  ph0.tmp = tmp1/sum(tmp1)

  ph0 = ph0.tmp[1]  # prob of the null
  pp0dash = ph0.tmp[-1]  # pps of variants

  # posterior probs of true system
  pp.tmp = exp(bf - coloc:::logsum(bf))
  pp0 = pp.tmp/sum(pp.tmp)

  # estimate mu
  muhat = mean(c(sum(abs(z) * pp0dash), (1 - ph0) * max(abs(z))))

  #### corrected coverage
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal
  # simulate ERR matrix

  ERR = mvtnorm:::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
  pp_ERR = function(Zj) {
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    zstar = mexp.zm + ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = coloc:::logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp/rowSums(pp.tmp)
  }
  # simulate pp systems
  pps <- mapply(pp_ERR, zj, SIMPLIFY = FALSE)
  # consider different CV as causal in each list
  n_pps <- length(pps)
  args <- 1:nsnps

  # obtain credible set for each simulation
  d5 <- lapply(1:n_pps, function(x) {
    credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
  })

  d5_trim <- lapply(d5, function(p) p[which(p$nvar==nvar),])

  propor_cov <- lapply(d5_trim, prop_cov) %>% unlist()

  nsims <- lapply(d5_trim, function(x) dim(x)[1]) %>% unlist()

  contained <- lapply(d5_trim, function(p) p$covered) %>% unlist()

  pp.vec <- rep(pp0, times=nsims)

  sum(contained * pp.vec)/sum(pp.vec)
}

#' Obtain corrected coverage estimate using estimated effect sizes and their standard errors (limiting simulations used for estimation to those with correct nvar)
#'
#' This function requires the marginal summary statistics from GWAS and an nvar value. It should only be used when nvar is very low ($<3$) and there is some evidence to suggest that only simulated credible sets with this nvar value should be used to derive the corrected coverage estimate.
#' @rdname corrcov_nvar_bhat
#' @title Obtain corrected coverage estimate using estimated effect sizes and their standard errors
#' @param bhat Estimated effect sizes from single-SNP logistic regressions
#' @param V Variance of estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param nvar The number of variants that simulated credible sets used for estimation should contain
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 10000 default because we trim out the ones without correct nvar)
#' @return Corrected coverage estimate
#' @export
#'
corrcov_nvar_bhat <- function(bhat, V, N0, N1, Sigma, nvar, thr = 0.95, W = 0.2, nrep = 10000) {
  z <- bhat/sqrt(V)
  r <- W^2/(W^2 + V)
  bf = 0.5 * (log(1 - r) + (r * z^2))
  p1 = 1e-04  # hard code
  nsnps = length(bf)
  prior = c(1 - nsnps * p1, rep(p1, nsnps))
  tmp = c(1, bf)  # add on extra for null model
  my.denom = coloc:::logsum(tmp + prior)
  tmp1 = exp(tmp + prior - my.denom)
  ph0.tmp = tmp1/sum(tmp1)

  ph0 = ph0.tmp[1]  # prob of the null
  pp0dash = ph0.tmp[-1]  # pps of variants

  # posterior probs of true system
  pp.tmp = exp(bf - coloc:::logsum(bf))
  pp0 = pp.tmp/sum(pp.tmp)

  # estimate mu
  muhat = mean(c(sum(abs(z) * pp0dash), (1 - ph0) * max(abs(z))))

  #### corrected coverage
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal
  # simulate ERR matrix

  ERR = mvtnorm:::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
  pp_ERR = function(Zj) {
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    zstar = mexp.zm + ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = coloc:::logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp/rowSums(pp.tmp)
  }

  # simulate pp systems
  pps <- mapply(pp_ERR, zj, SIMPLIFY = FALSE)
  # consider different CV as causal in each list
  n_pps <- length(pps)
  args <- 1:nsnps

  # obtain credible set for each simulation
  d5 <- lapply(1:n_pps, function(x) {
    credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
  })

  d5_trim <- lapply(d5, function(p) p[which(p$nvar==nvar),])

  propor_cov <- lapply(d5_trim, prop_cov) %>% unlist()

  nsims <- lapply(d5_trim, function(x) dim(x)[1]) %>% unlist()

  contained <- lapply(d5_trim, function(p) p$covered) %>% unlist()

  pp.vec <- rep(pp0, times=nsims)

  sum(contained * pp.vec)/sum(pp.vec)
}

#' Obtain confidence interval for corrected coverage estimate using Z-scores and mafs
#'
#' @rdname corrcov_CI
#' @title Obtain confidence interval for corrected coverage estimate using Z-scores and mafs
#' @param z Marginal Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 1000 default)
#' @param CI The size of the confidence interval (as a decimal)
#' @export
#' @return Corrected coverage estimate
#'
corrcov_CI <- function(z, f, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000, CI = 0.95){
  corrcov_reps = replicate(100, corrcov(z, f, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000))
  quantile(corrcov_reps, probs = c((1-CI)/2, (CI+1)/2))
}

#' Obtain confidence interval for corrected coverage estimate using estimated effect sizes and their standard errors
#'
#' @rdname corrcov_CI_bhat
#' @title Obtain confidence interval for corrected coverage estimate using estimated effect sizes and their standard errors
#' @param bhat Estimated effect sizes from single-SNP logistic regressions
#' @param V Variance of estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 1000 default)
#' @param CI The size of the confidence interval (as a decimal)
#' @return Corrected coverage estimate
#' @export
#'
corrcov_CI_bhat <- function(bhat, V, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000, CI = 0.95){
  corrcov_reps = replicate(100, corrcov_bhat(bhat, V, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000))
  quantile(corrcov_reps, probs = c((1-CI)/2, (CI+1)/2))
}

