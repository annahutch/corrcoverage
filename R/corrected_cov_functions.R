#' Get credible set of variants
#'
#' If the CV parameter is supplied (index of causal variant) then the
#' output includes a binary indicator of whether the CV is contained in the set
#' @title Get credible set of variants
#' @param pp Vector of posterior probabilities of causality
#' @param CV Optional parameter: Index of CV
#' @param thr Minimum threshold for credible set size (default is 0.95)
#' @export
#' @return list of the variants in the credible set, the claimed.cov (cumulative sum of the posterior probabilities of the variants forming the credible set), binary covered indicator (1 if CV is contained in the credible set) and nvar (number of variants in the set)
credset <- function(pp, CV, thr = 0.95) {
    o <- order(pp, decreasing = TRUE)  # order index for true pp
    cumpp <- cumsum(pp[o])  # cum sums of ordered pps
    wh <- which(cumpp > thr)[1]  # how many needed to exceed thr
    size <- cumpp[wh]
    if (missing(CV)) {
        data.frame(claimed.cov = size, nvar = wh)
    } else {
        contained = as.numeric(CV %in% o[1:wh])
        list(credset = o[1:wh], claimed.cov = size, covered = contained, nvar = wh)
    }
}

#' Quicker credset function for matrix of posterior probabilities (using RCpp)
#'
#' @title Get credible set of variants from matrix of pps (Rcpp)
#' @param pp Matrix of posterior probabilities of causality (one row per system)
#' @param CV Vector of CV indices (one per system/row)
#' @param thr Minimum threshold for credible set size (default is 0.95)
#'
#' @return Data.frame of claimed coverage (sum of posterior probabilities of variants in the set), binary covered indicator and number of variants (nvar).
#' @useDynLib corrcoverage
#' @importFrom Rcpp sourceCpp
#' @export
credsetC <- function(pp, CV = iCV, thr = 0.95) {
    ret <- credsetmat(pp, CV, thr)  ## list 1 = wh, 2 = size, 3=contained
    data.frame(claimed.cov = ret[[2]], covered = ret[[3]], nvar = ret[[1]])
}

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

