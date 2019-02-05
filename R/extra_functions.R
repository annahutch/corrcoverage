#' Variance of the estimated effect size for case-control
#'
#' @title Var.data.cc
#' @param f Minor allele frequencies
#' @param N Total sample size
#' @param s Proportion of cases (N1/N0+N1)
#' @return Variance of estimated effect size, V.
#' @author Claudia Giambartolomei
Var.data.cc <- function(f, N, s) {
    1/(2 * N * f * (1 - f) * s * (1 - s))
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
    my.res <- my.max + log(sum(exp(x - my.max)))
    return(my.res)
}

#' Correlation matrix of SNPS
#'
#' @title cor2
#' @param x Phased haplotype matrix, rows as samples and columns as SNPs
#' @return Correlation matrix
#' @author Chris Wallace
cor2 <- function(x) {
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
    file.ldd = "/home/cew54/share/Data/reference/lddetect/EUR/fourier_ls-chr22.bed"
    file.vcf = "/home/cew54/share/Data/reference/UK10K/BCF/chr22.bcf.gz"
    
    ## ldblocks
    ldd <- fread(file.ldd)
    
    ## split bcf by ldblocks
    ldd[, `:=`(blocknum, 1:.N)]
    ldd[, `:=`(dist, stop - start)]
    ldd[, `:=`(comm, paste0("/home/cew54/share/bin/bcftools view ", file.vcf, " --min-af 0.02:minor --max-alleles 2 --min-alleles 2 ", " -r chr", 
        22, ":", start, "-", stop, " -Ov "))]  # -o ',tmp)]
    gethap <- function(i) {
        y = fread(ldd$comm[i])
        ha <- simGWAS:::vcf2haps(as.matrix(y[, -c(1:9)]))
        rownames(ha) <- paste0("pos", y$POS)
        t(ha)
    }
    
    block <- sample(which(ldd$dist < 1200000), 1)  # use smallest LD block to be fast
    h <- gethap(block)  # rows=samples, cols=snps
    use <- apply(h, 2, var) > 0 & colMeans(h) > 0.01 & colMeans(h) < 0.99  # no monomorphs
    h <- h[, use, drop = FALSE]
    nmax <- floor(nsnp/2)  # keep this small to make simulations fast
    if (ncol(h) > nmax) {
        start <- sample(1:ncol(h), 1)  # choose a random starting point
        start2 <- start + floor(ncol(h)/2)  # choose another random starting point, the other end
        idx <- c(start + 1:nmax, start2 + 1:nmax)  # doing this means we're not just looking at one LD block of SNPs, so there will be some snps completely seperate to CV
        idx <- idx%%ncol(h)
        h <- h[, idx]
    }
    list(h = h, block = block)
}

#' Build and predict corrected coverage using logistic GAM
#'
#' @rdname pred_logit
#' @title pred_logit
#' @param x data.frame with column for 'covered' and 'logit.claim'
#' @return Predicted probability of covered
pred_logit <- function(x) {
    m <- mgcv:::gam(covered ~ s(logit.claim), data = x, family = "binomial")
    invlogit(predict(m, newdata = data.frame(logit.claim = logit(claim0))))
}

#' Prediction for corrected coverage if gam cannot be fitted
#'
#' @rdname pred_na
#' @title pred_na
#' @param x data.frame with a binary 'covered' column
#' @return Predicted coverage
#' @author Anna Hutchinson
pred_na <- function(x) {
    mean(x$covered)
}

#' Estimate for true effect at CV
#'
#' @rdname mu.est
#' @title mu_est
#' @return Estimate of true effect at CV
mu_est <- function(X) {
    y = seq(0, 20, 0.005)
    x <- sapply(y, function(m) mean(abs(rnorm(50000, mean = m))))
    approx(x, y, xout = X)$y
}
