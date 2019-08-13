#' Corrected coverage estimate of the causal variant in the credible set
#'
#' Requires an estimate of the true effect at the CV (e.g. use maximum absolute z-score or output from corrcoverage::mu_est function)
#' @rdname corrected_cov
#' @title Corrected coverage estimate of the causal variant in the credible set
#' @param pp0 Posterior probabilities of SNPs
#' @param mu The true effect at the CV (estimate using corrcoverage::mu_est function)
#' @param V Variance of the estimated effect size (can be obtained using coloc::Var.beta.cc function)
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (0.95 default)
#' @param W Prior for the standard deviation of the effect size parameter, beta (W=0.2 default)
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
#' maf <- colMeans(X)
#'
#' ## generate V (variance of estimated effect sizes)
#' varbeta <- Var.data.cc(f = maf, N = 5000, s = 0.5)
#'
#' pp <- rnorm(nsnps, 0.2, 0.05)
#' pp <- pp/sum(pp)
#'
#' corrected_cov(pp0 = pp, mu = 4, V = varbeta, Sigma = LD, thr = 0.95, nrep = 100)
#'
#' @export
#' @author Anna Hutchinson
corrected_cov <- function(pp0, mu, V, Sigma, thr = 0.95, W = 0.2, nrep = 1000) {

    nsnps = length(pp0)
    temp = diag(x = mu, nrow = nsnps, ncol = nsnps)
    zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal

    # simulate ERR matrix
    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)

    # calculate r
    r = W^2/(W^2 + V)

    # simulate pp systems
    pps = mapply(.zj_pp, Zj = zj, MoreArgs = list(int.Sigma = Sigma, int.nrep = nrep, int.ERR = ERR, int.r = r), SIMPLIFY =     FALSE)

    # consider different CV as causal in each list
    n_pps <- length(pps)
    args <- 1:nsnps

    # obtain credible set for each simulation
    d5 <- lapply(1:n_pps, function(x) {
        credsetC(pps[[x]], CV = rep(args[x], dim(pps[[x]])[1]), thr = thr)
    })

    propcov <- lapply(d5, prop_cov) %>% unlist()
    sum(propcov * pp0)
}

#' Corrected coverage estimate of the causal variant in the credible set
#'
#' Requires an estimate of the true effect at the CV (e.g. use maximum absolute z-score or output from corrcoverage::mu_est function)
#' @rdname corrected_cov
#' @title Corrected coverage estimate of the causal variant in the credible set
#' @param pp0 Posterior probabilities of SNPs
#' @param mu The true effect at the CV (estimate using corrcoverage::mu_est function)
#' @param V Variance of the estimated effect size (can be obtained using coloc::Var.beta.cc function)
#' @param Sigma SNP correlation matrix
#' @param thr Minimum threshold for fine-mapping experiment (0.95 default)
#' @param W Prior for the standard deviation of the effect size parameter, beta (W=0.2 default)
#' @param nrep Number of posterior probability systems to simulate for each variant considered causal (nrep = 1000 default)
#' @param pp0min only average over SNPs with pp0 > pp0min
#' @return Corrected coverage estimate
#'
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#' N0 <- 5000
#' N1 <- 5000
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
#' maf <- colMeans(X)
#'
#' ## generate V (variance of estimated effect sizes)
#' varbeta <- Var.data.cc(f = maf, N = 5000, s = 0.5)
#'
#' pp <- rnorm(nsnps, 0.2, 0.05)
#' pp <- pp/sum(pp)
#'
#' corrected_cov(pp0 = pp, mu = 4, V = varbeta, Sigma = LD, thr = 0.95, nrep = 100)
#'
#' @export
#' @author Anna Hutchinson
corrected_cov_faster <- function(pp0, mu, V, Sigma, thr = 0.95, W = 0.2, nrep = 1000, pp0min=0.001) {

    nsnps = length(pp0)
    temp = diag(x = mu, nrow = nsnps, ncol = nsnps)
    usesnps=which(pp0 > pp0min)
    zj = lapply(usesnps, function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal
    
    # simulate ERR matrix
    ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)

    # calculate r
    r = W^2/(W^2 + V)

    # simulate pp systems
    pps = mapply(.zj_pp, Zj = zj, MoreArgs = list(int.Sigma = Sigma, int.nrep = nrep, int.ERR = ERR, int.r = r), SIMPLIFY =     FALSE)

    # consider different CV as causal in each list
    n_pps <- length(pps)
    args <- 1:nsnps

    # obtain credible set for each simulation
    d5 <- lapply(1:n_pps, function(x) {
        credsetC(pps[[x]], CV = rep(usesnps[x], dim(pps[[x]])[1]), thr = thr)
    })

    propcov <- lapply(d5, prop_cov) %>% unlist()
    sum(propcov * pp0[usesnps])
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
#' @param thr Minimum threshold for fine-mapping experiment (default 0.95)
#' @param W Prior for the standard deviation of the effect size parameter, beta (default 0.2)
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (default 1000)
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
#' ## generate example LD matrix
#' library(mvtnorm)
#' nsamples = 1000
#'
#' simx <- function(nsnps, nsamples, S, maf=0.1) {
#'     mu <- rep(0,nsnps)
#'     rawvars <- rmvnorm(n=nsamples, mean=mu, sigma=S)
#'     pvars <- pnorm(rawvars)
#'     x <- qbinom(1-pvars, 1, maf)
#' }
#'
#' S <- (1 - (abs(outer(1:nsnps,1:nsnps,`-`))/nsnps))^4
#' X <- simx(nsnps,nsamples,S)
#' LD <- cor2(X)
#' maf <- colMeans(X)
#'
#' corrcov(z = z_scores, f = maf, N0, N1, Sigma = LD, thr = 0.95)
#'
#' @export
#' @author Anna Hutchinson
corrcov <- function(z, f, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000) {

    varbeta = 1/(2 * (N0 + N1) * f * (1 - f) * (N1/(N0 + N1)) * (1 - (N1/(N0 + N1))))

    pp = ppfunc(z, V = varbeta, W = 0.2)

    muhat = sum(abs(z) * pp)

    corrected_cov(pp0 = pp, mu = muhat, V = varbeta, Sigma, thr, W, nrep)
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
#' @param thr Minimum threshold for fine-mapping experiment (default 0.95)
#' @param W Prior for the standard deviation of the effect size parameter, beta (default 0.2)
#' @param nrep The number of simulated posterior probability systems to consider for the corrected coverage estimate (default 1000)
#' @return Corrected coverage estimate
#'
#' @examples
#'
#' set.seed(1)
#' nsnps <- 100
#' N0 <- 1000 # number of controls
#' N1 <- 1000 # number of cases
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
#' maf <- colMeans(X)
#'
#' varbeta <- Var.data.cc(f = maf, N = N0 + N1, s = N1/(N0+N1))
#'
#' bhats = rnorm(nsnps, 0, 0.2) # log OR
#'
#' corrcov_bhat(bhat = bhats, V = varbeta, N0, N1, Sigma = LD)
#'
#' @export
#' @author Anna Hutchinson
corrcov_bhat <- function(bhat, V, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000) {

    z = bhat/sqrt(V)

    pp = ppfunc(z, V, W = 0.2)

    muhat = sum(abs(z) * pp)

    corrected_cov(pp0 = pp, mu = muhat, V, Sigma, thr, W, nrep)
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
#' @param thr Minimum threshold for fine-mapping experiment (default 0.95)
#' @param W Prior for the standard deviation of the effect size parameter, beta (default 0.2)
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
#' maf <- colMeans(X)
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

  pp = ppfunc(z, V = varbeta, W = 0.2)

  muhat = sum(abs(z) * pp)

  nsnps = length(pp)

  #### corrected coverage

  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal

  # simulate ERR matrix

  ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)

  r = W^2/(W^2 + varbeta)

  pps = mapply(.zj_pp, Zj = zj, MoreArgs = list(int.Sigma = Sigma, int.nrep = nrep, int.ERR = ERR, int.r = r), SIMPLIFY =     FALSE)

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

  pp.vec <- rep(pp, times=nsims)

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
#' @param thr Minimum threshold for fine-mapping experiment (default 0.95)
#' @param W Prior for the standard deviation of the effect size parameter, beta (default 0.2)
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
#' maf <- colMeans(X)
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

  pp = ppfunc(z, V, W = 0.2)

  muhat = sum(abs(z) * pp)

  nsnps = length(pp)

  #### corrected coverage

  temp = diag(x = muhat, nrow = nsnps, ncol = nsnps)
  zj = lapply(seq_len(nrow(temp)), function(i) temp[i, ])  # nsnp zj vectors for each snp considered causal

  # simulate ERR matrix

  ERR = mvtnorm::rmvnorm(nrep, rep(0, ncol(Sigma)), Sigma)

  r = W^2/(W^2 + V)

  pps = mapply(.zj_pp, Zj = zj, MoreArgs = list(int.Sigma = Sigma, int.nrep = nrep, int.ERR = ERR, int.r = r), SIMPLIFY =     FALSE)

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

  pp.vec <- rep(pp, times=nsims)

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
#' @param thr Minimum threshold for fine-mapping experiment (default 0.95)
#' @param W Prior for the standard deviation of the effect size parameter, beta (default 0.2)
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
#' ## generate example LD matrix
#' library(mvtnorm)
#' nsamples = 1000
#'
#' simx <- function(nsnps, nsamples, S, maf=0.1) {
#'     mu <- rep(0,nsnps)
#'     rawvars <- rmvnorm(n=nsamples, mean=mu, sigma=S)
#'     pvars <- pnorm(rawvars)
#'     x <- qbinom(1-pvars, 1, maf)
#' }
#'
#' S <- (1 - (abs(outer(1:nsnps,1:nsnps,`-`))/nsnps))^4
#' X <- simx(nsnps,nsamples,S)
#' LD <- cor2(X)
#'
#' corrcov_CI(z = z_scores, f = maf, N0, N1, Sigma = LD)
#' }
#'
#' @export
#'
#' @author Anna Hutchinson
corrcov_CI <- function(z, f, N0, N1, Sigma, thr = 0.95, W = 0.2, nrep = 1000, CI = 0.95){
  corrcov_reps = replicate(100, corrcov(z, f, N0, N1, Sigma, thr, W, nrep))
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
  corrcov_reps = replicate(100, corrcov_bhat(bhat, V, N0, N1, Sigma, thr, W, nrep))
  stats::quantile(corrcov_reps, probs = c((1-CI)/2, (CI+1)/2))
}
