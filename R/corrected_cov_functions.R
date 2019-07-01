#' Corrected coverage estimate of the causal variant in the credible set
#'
#' Requires an estimate of the true effect at the CV (e.g. use maximum absolute z-score or output from corrcoverage::mu_est function)
#' @rdname corrected_cov
#' @title Corrected coverage estimate of the causal variant in the credible set
#' @param mu The true effect at the CV (estimate using corrcoverage::mu_est function)
#' @param V Variance of the estimated effect size (can be obtained using Var.beta.cc function)
#' @param W Prior for the standard deviation of the effect size parameter, beta (W=0.2 default)
#' @param Sigma SNP correlation matrix
#' @param pp0 Posterior probabilities of SNPs
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param nrep Number of posterior probability systems to simulate for each variant considered causal (nrep = 1000 default)
#' @return Corrected coverage estimate
#'
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#' N0 <- 5000
#' N1 <- 5000
#'
#' ## generate example LD matrix (https://chr1swallace.github.io/simGWAS/articles/intro.html)
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' LD <- cor2(haps)
#' maf <- colMeans(haps)
#'
#' ## generate V (variance of estimated effect sizes)
#' varbeta <- Var.data.cc(f = maf, N = 5000, s = 0.5)
#'
#' pp <- rnorm(nsnps, 0.2, 0.05)
#' pp <- pp/sum(pp)
#'
#' corrected_cov(mu = 4, V = varbeta, Sigma = LD, pp0 = pp, thr = 0.95, nrep = 100)
#'
#' @export
#' @author Anna Hutchinson
corrected_cov <- function(mu, V, W = 0.2, Sigma, pp0, thr = 0.95, nrep = 1000) {

    nsnps = length(pp0)
    temp = diag(x = mu, nrow = nsnps, ncol = nsnps)
    zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal

    # simulate ERR matrix
    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    r = W^2/(W^2 + V)

    pp_ERR = function(Zj) {
        exp.zm = Zj %*% Sigma
        mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
        zstar = mexp.zm + ERR
        bf = 0.5 * (log(1 - r) + (r * zstar^2))
        denom = logsum(bf)
        pp.tmp = exp(bf - denom)  # convert back from log scale
        pp.tmp/rowSums(pp.tmp)
    }

    # simulate pp systems
    pps = mapply(pp_ERR, zj, SIMPLIFY = FALSE)

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


#' Corrected coverage estimate using Z-scores and mafs
#'
#' This function only requires the marginal summary statistics from GWAS
#' @rdname corrcov
#' @title Corrected coverage estimate using Z-scores and mafs
#' @param z Marginal Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 1000 default)
#' @return Corrected coverage estimate
#'
#' @examples
#'
#' set.seed(1)
#' nsnps = 100
#' N0 = 5000
#' N1 = 5000
#' z_scores <- rnorm(nsnps, 0, 3) # simulate a vector of Z-scores
#'
#' # simulate fake haplotypes to obtain MAFs and LD matrix
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' maf <- colMeans(haps)
#' LD <- cor2(haps)
#'
#' corrcov(z = z_scores, f = maf, N0, N1, Sigma = LD, thr = 0.95)
#'
#' @export
#' @author Anna Hutchinson
corrcov <- function(z, f, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000) {
    varbeta = 1/(2 * (N0 + N1) * f * (1 - f) * (N1/(N0 + N1)) * (1 - (N1/(N0 + N1))))
    r = W^2/(W^2 + varbeta)
    bf = 0.5 * (log(1 - r) + (r * z^2))
    p1 = 1e-04
    nsnps = length(bf)
    prior = c(1 - nsnps * p1, rep(p1, nsnps))
    tmp = c(1, bf)  # add on extra for null model
    my.denom = logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    ph0.tmp = tmp1/sum(tmp1)

    ph0 = ph0.tmp[1]  # prob of the null
    pp0dash = ph0.tmp[-1]  # pps of variants

    # posterior probs of true system
    pp.tmp = exp(bf - logsum(bf))
    pp0 = pp.tmp/sum(pp.tmp)

    # estimate mu
    muhat = mean(c(sum(abs(z) * pp0dash), (1 - ph0) * max(abs(z))))

    #### corrected coverage
    temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
    zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal
    # simulate ERR matrix

    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    pp_ERR = function(Zj) {
        exp.zm = Zj %*% Sigma
        mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
        zstar = mexp.zm + ERR
        bf = 0.5 * (log(1 - r) + (r * zstar^2))
        denom = logsum(bf)
        pp.tmp = exp(bf - denom)  # convert back from log scale
        pp.tmp/rowSums(pp.tmp)
    }
    # simulate pp systems
    pps = mapply(pp_ERR, zj, SIMPLIFY = FALSE)
    # consider different CV as causal in each list
    n_pps = length(pps)
    args = 1:nsnps

    # obtain credible set for each simulation
    d5 <- lapply(1:n_pps, function(x) {
        credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })

    prop_cov = lapply(d5, prop_cov) %>% unlist()

    sum(prop_cov * pp0)
}

#' Corrected coverage estimate using estimated effect sizes and their standard errors
#'
#' This function only requires the marginal summary statistics from GWAS
#' @rdname corrcov_bhat
#' @title Corrected coverage estimate using estimated effect sizes and their standard errors
#' @param bhat Estimated effect sizes from single-SNP logistic regressions
#' @param V Variance of estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 1000 default)
#' @return Corrected coverage estimate
#'
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#' N0 <- 1000 # number of controls
#' N1 <- 1000 # number of cases
#'
#' # simulate fake haplotypes to obtain MAFs and LD matrix
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' maf <- colMeans(haps)
#' LD <- cor2(haps)
#'
#' varbeta <- Var.data.cc(f = maf, N = N0 + N1, s = N1/(N0+N1))
#'
#' bhats = rnorm(nsnps, 0, 0.2) # log OR
#'
#'
#' corrcov_bhat(bhat = bhats, V = varbeta, N0, N1, Sigma = LD)
#'
#' @export
#' @author Anna Hutchinson
corrcov_bhat <- function(bhat, V, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000) {
    z = bhat/sqrt(V)
    r = W^2/(W^2 + V)
    bf = 0.5 * (log(1 - r) + (r * z^2))
    p1 = 1e-04
    nsnps = length(bf)
    prior = c(1 - nsnps * p1, rep(p1, nsnps))
    tmp = c(1, bf)  # add on extra for null model
    my.denom = logsum(tmp + prior)
    tmp1 = exp(tmp + prior - my.denom)
    ph0.tmp = tmp1/sum(tmp1)

    ph0 = ph0.tmp[1]  # prob of the null
    pp0dash = ph0.tmp[-1]  # pps of variants

    # posterior probs of true system
    pp.tmp = exp(bf - logsum(bf))
    pp0 = pp.tmp/sum(pp.tmp)

    # estimate mu
    muhat = mean(c(sum(abs(z) * pp0dash), (1 - ph0) * max(abs(z))))

    #### corrected coverage
    temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
    zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal

    # simulate ERR matrix

    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
    pp_ERR = function(Zj) {
        exp.zm = Zj %*% Sigma
        mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
        zstar = mexp.zm + ERR
        bf = 0.5 * (log(1 - r) + (r * zstar^2))
        denom = logsum(bf)
        pp.tmp = exp(bf - denom)  # convert back from log scale
        pp.tmp/rowSums(pp.tmp)
    }

    # simulate pp systems
    pps = mapply(pp_ERR, zj, SIMPLIFY = FALSE)
    # consider different CV as causal in each list
    n_pps = length(pps)
    args = 1:nsnps

    # obtain credible set for each simulation
    d5 <- lapply(1:n_pps, function(x) {
        credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })

    prop_cov = lapply(d5, prop_cov) %>% unlist()

    sum(prop_cov * pp0)
}

#' Obtain corrected coverage estimate using Z-scores and mafs (limiting simulations used for estimation to those with correct nvar)
#'
#' This function requires the marginal summary statistics from GWAS and an nvar value. It should only be used when nvar is very low ($<3$) and there is some evidence to suggest that only simulated credible sets with this nvar value should be used to derive the corrected coverage estimate.
#' @rdname corrcov_nvar
#' @title Corrected coverage estimate using Z-scores and mafs (fixing nvar)
#' @param z Marginal Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param nvar The number of variants that simulated credible sets used for estimation should contain
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected
#'            coverage estimate (nrep = 10000 default because we trim out the ones without correct
#'            nvar so need this to be hugh)
#' @return Corrected coverage estimate
#'
#' @examples
#'
#' set.seed(1)
#' nsnps = 100
#' N0 = 5000
#' N1 = 5000
#' z_scores <- rnorm(nsnps, 0, 3) # simulate a vector of Z-scores
#'
#' # simulate fake haplotypes to obtain MAFs and LD matrix
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' maf <- colMeans(haps)
#' LD <- cor2(haps)
#'
#' corrcov_nvar(z = z_scores, f = maf, N0, N1, Sigma = LD, nvar = 1, nrep = 100)
#'
#' # note that nrep should be at least the default value (nrep = 10000) but is
#' # lower here for speed of computation
#'
#' @export

#' @author Anna Hutchinson
corrcov_nvar <- function(z, f, N0, N1, Sigma, nvar, thr = 0.95, W = 0.2, nrep = 10000) {
  varbeta = 1/(2 * (N0 + N1) * f * (1 - f) * (N1/(N0 + N1)) * (1 - (N1/(N0 + N1))))
  r = W^2/(W^2 + varbeta)
  bf = 0.5 * (log(1 - r) + (r * z^2))
  p1 = 1e-04  # hard code
  nsnps = length(bf)
  prior = c(1 - nsnps * p1, rep(p1, nsnps))
  tmp = c(1, bf)  # add on extra for null model
  my.denom = logsum(tmp + prior)
  tmp1 = exp(tmp + prior - my.denom)
  ph0.tmp = tmp1/sum(tmp1)

  ph0 = ph0.tmp[1]  # prob of the null
  pp0dash = ph0.tmp[-1]  # pps of variants

  # posterior probs of true system
  pp.tmp = exp(bf - logsum(bf))
  pp0 = pp.tmp/sum(pp.tmp)

  # estimate mu
  muhat = mean(c(sum(abs(z) * pp0dash), (1 - ph0) * max(abs(z))))

  #### corrected coverage
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal
  # simulate ERR matrix

  ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
  pp_ERR = function(Zj) {
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    zstar = mexp.zm + ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
    pp.tmp = exp(bf - denom)  # convert back from log scale
    pp.tmp/rowSums(pp.tmp)
  }
  # simulate pp systems
  pps = mapply(pp_ERR, zj, SIMPLIFY = FALSE)
  # consider different CV as causal in each list
  n_pps = length(pps)
  args = 1:nsnps

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
#' @title Corrected coverage estimate using estimated effect sizes and their standard errors (fixing nvar)
#' @param bhat Estimated effect sizes from single-SNP logistic regressions
#' @param V Variance of estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param nvar The number of variants that simulated credible sets used for estimation should contain
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected
#'            coverage estimate (nrep = 10000 default because we trim out the ones without correct
#'            nvar so need this to be hugh)
#' @return Corrected coverage estimate
#'
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#' N0 <- 5000 # number of controls
#' N1 <- 5000 # number of cases
#'
#' # simulate fake haplotypes to obtain MAFs and LD matrix
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' maf <- colMeans(haps)
#' LD <- cor2(haps)
#'
#' varbeta <- Var.data.cc(f = maf, N = N0 + N1, s = N1/(N0+N1))
#'
#' bhats = rnorm(nsnps,0,0.2) # log OR
#'
#' corrcov_nvar_bhat(bhat = bhats, V = varbeta, N0, N1, Sigma = LD, nvar = 1, nrep = 1000)
#'
#' # note that nrep should be at least the default value (nrep = 10000) but is
#' # lower here for speed of computation
#'
#' @export
#'
#' @author Anna Hutchinson
corrcov_nvar_bhat <- function(bhat, V, N0, N1, Sigma, nvar, thr = 0.95, W = 0.2, nrep = 10000) {
  z = bhat/sqrt(V)
  r = W^2/(W^2 + V)
  bf = 0.5 * (log(1 - r) + (r * z^2))
  p1 = 1e-04  # hard code
  nsnps = length(bf)
  prior = c(1 - nsnps * p1, rep(p1, nsnps))
  tmp = c(1, bf)  # add on extra for null model
  my.denom = logsum(tmp + prior)
  tmp1 = exp(tmp + prior - my.denom)
  ph0.tmp = tmp1/sum(tmp1)

  ph0 = ph0.tmp[1]  # prob of the null
  pp0dash = ph0.tmp[-1]  # pps of variants

  # posterior probs of true system
  pp.tmp = exp(bf - logsum(bf))
  pp0 = pp.tmp/sum(pp.tmp)

  # estimate mu
  muhat = mean(c(sum(abs(z) * pp0dash), (1 - ph0) * max(abs(z))))

  #### corrected coverage
  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal
  # simulate ERR matrix

  ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)
  pp_ERR = function(Zj) {
    exp.zm = Zj %*% Sigma
    mexp.zm = matrix(exp.zm, nrep, length(Zj), byrow = TRUE)  # matrix of Zj replicated in each row
    zstar = mexp.zm + ERR
    bf = 0.5 * (log(1 - r) + (r * zstar^2))
    denom = logsum(bf)  # logsum(x) = max(x) + log(sum(exp(x - max(x)))) so sum is not inf
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
#' @title Confidence interval for corrected coverage estimate using Z-scores and mafs
#' @param z Marginal Z-scores
#' @param f Minor allele frequencies
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 1000 default)
#' @param CI The size of the confidence interval (as a decimal)
#' @return CI for corrected coverage estimate
#'
#' @examples
#'
#' \dontrun{
#'
#'  # this is a long running example
#' set.seed(1)
#' nsnps = 100
#' N0 = 5000
#' N1 = 5000
#' z_scores <- rnorm(nsnps, 0, 3) # simulate a vector of Z-scores
#'
#' # simulate fake haplotypes to obtain MAFs and LD matrix
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' maf <- colMeans(haps)
#' LD <- cor2(haps)
#'
#' corrcov_CI(z = z_scores, f = maf, N0, N1, Sigma = LD)
#' }
#'
#' @export
#'
#' @author Anna Hutchinson
corrcov_CI <- function(z, f, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000, CI = 0.95){
  corrcov_reps = replicate(100, corrcov(z, f, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000))
  stats::quantile(corrcov_reps, probs = c((1-CI)/2, (CI+1)/2))
}

#' Obtain confidence interval for corrected coverage estimate using estimated effect sizes and their standard errors
#'
#' @rdname corrcov_CI_bhat
#' @title Confidence interval for corrected coverage estimate using estimated effect sizes and their standard errors
#' @param bhat Estimated effect sizes from single-SNP logistic regressions
#' @param V Variance of estimated effect sizes
#' @param N0 Number of controls
#' @param N1 Number of cases
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (default is 0.95)
#' @param W Prior for the standard deviation of the effect size parameter beta
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (nrep = 1000 default)
#' @param CI The size of the confidence interval (as a decimal)
#' @return CI for corrected coverage estimate
#'
#' @examples
#'
#' \dontrun{
#'  # this is a long running example
#' set.seed(1)
#' nsnps <- 100
#' N0 <- 5000 # number of controls
#' N1 <- 5000 # number of cases
#'
#' # simulate fake haplotypes to obtain MAFs and LD matrix
#' nhaps <- 1000
#' lag <- 5
#' maf.tmp <- runif(nsnps+lag, 0.05, 0.5) # common SNPs
#' laghaps <- do.call("cbind", lapply(maf.tmp, function(f) rbinom(nhaps,1,f)))
#' haps <- laghaps[,1:nsnps]
#' for(j in 1:lag)
#'    haps <- haps + laghaps[,(1:nsnps)+j]
#' haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
#' maf <- colMeans(haps)
#' LD <- cor2(haps)
#'
#' varbeta <- Var.data.cc(f = maf, N = N0 + N1, s = N1/(N0+N1))
#'
#' bhats = rnorm(nsnps,0,0.2) # log OR
#'
#' corrcov_CI_bhat(bhat = bhats, V = varbeta, N0, N1, Sigma = LD)
#' }
#'
#' @export
#'
#' @author Anna Hutchinson
corrcov_CI_bhat <- function(bhat, V, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000, CI = 0.95){
  corrcov_reps = replicate(100, corrcov_bhat(bhat, V, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000))
  stats::quantile(corrcov_reps, probs = c((1-CI)/2, (CI+1)/2))
}

