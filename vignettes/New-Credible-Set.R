## ----setup, set.seed = 18, include=FALSE---------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(18)
library(corrcoverage)

## ------------------------------------------------------------------------
library(simGWAS)

#  Simulate reference haplotypes
nsnps <- 200
nhaps <- 1000
lag <- 30  # genotypes are correlated between neighbouring variants
maf <- runif(nsnps + lag, 0.05, 0.5)  # common SNPs
laghaps <- do.call("cbind", lapply(maf, function(f) rbinom(nhaps, 1, f)))
haps <- laghaps[, 1:nsnps]
for (j in 1:lag) haps <- haps + laghaps[, (1:nsnps) + j]
haps <- round(haps/matrix(apply(haps, 2, max), nhaps, nsnps, byrow = TRUE))
snps <- colnames(haps) <- paste0("s", 1:nsnps)
freq <- as.data.frame(haps + 1)
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)
MAF <- colMeans(freq[, snps] - 1)  # minor allele frequencies
CV <- sample(snps[which(colMeans(haps) > 0.1)], 1)
iCV <- sub("s", "", CV)  # index of cv
LD <- cor2(haps) # correlation between SNPs

## ------------------------------------------------------------------------
OR <- 1.1 # odds ratios
N0 <- 10000 # number of controls
N1 <- 10000 # number of cases
  
z0 <- simulated_z_score(N0 = N0, # number of controls
                        N1 = N1, # number of cases
                        snps = snps, # column names in freq of SNPs for which Z scores should be generated
                        W = CV, # causal variants, subset of snps
                        gamma.W = log(OR), # log odds ratios
                        freq = freq) # reference haplotypes

## ------------------------------------------------------------------------
varbeta <- Var.data.cc(f = MAF, N = N1+N0, s = N1/(N0+N1)) # variance of estimated effect size

## ------------------------------------------------------------------------
postprobs <- ppfunc(z = z0, V = varbeta) 

## ------------------------------------------------------------------------
muhat <- est_mu(z0, MAF, N0, N1)
muhat

