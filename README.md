
<!-- README.md is generated from README.Rmd. Please edit that file -->
corrcoverage
============

The `corrcoverage` R package can be used to obtain an accurate coverage estimate of the causal variant in a credible set that has been acquired using the Bayesian approach for fine-mapping ([Maller et al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/23104008), [Wakefield, 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359).).

Webpage: <https://annahutch.github.io/corrcoverage/>

------------------------------------------------------------------------

Installation
------------

You can install the released version of `corrcoverage` from [github](https://github.com/) with:

``` r
install.packages("devtools") # if not already installed
devtools::install_github("annahutch/corrcoverage")
```

------------------------------------------------------------------------

Examples
--------

For examples and usage, please see the vignettes \[<http://annahutch.github.io/corrcoverage/articles/my-vignette.html>\]

------------------------------------------------------------------------

Conversion Functions
--------------------

The `corrcoverage` package includes useful functions for converting between p-values, *Z*-scores, asymptotic Bayes factors and posterior probabilities of causality. The following table summarises these functions.

| Include null model? [1] | Input                     | Output                  | Function      |
|-------------------------|---------------------------|-------------------------|---------------|
| YES                     | p-values                  | log(ABF)                | `approx.bf.p` |
| YES                     | p-values                  | Posterior Probabilities | `pvals_pp`    |
| YES                     | Marginal *Z*-scores       | Posterior Probabilties  | `z0_pp`       |
| NO                      | Marginal *Z*-scores       | Posterior Probabilties  | `ppfunc`      |
| NO                      | Marginal *Z*-score matrix | Posterior Probabilties  | `z0func.mat`  |

Functions are also provided to simulate marginal *Z*-scores from joint *Z*-scores. The joint *Z*-scores are all 0, except at the causal variant where it is the 'true effect', *μ*. *μ* can be estimated using `est_mu` function which required sample sizes, *Z*-scores and minor allele frequencies.

`z_sim` simulates marginal *Z*-scores from joint *Z*-scores, whilst `zj_pp` can be used to simulate posterior probability systems from a joint *Z*-score vector. These functions work by first calculating *E*(*Z*<sub>*m*</sub>),
*E*(*Z*<sub>*m*</sub>)=*Z*<sub>*j*</sub> × *Σ*
 where *Σ* is the correlation matrix between SNPs.

We can then simulate more *Z*-score systems from a multivariate normal disribution with mean *E*(*Z*<sub>*m*</sub>) and variance *Σ*. This is a key step in our corrected coverage method.

------------------------------------------------------------------------

Our correction factor works by simulating many more credible sets from 'the same system as the original' and finding what proportion of these contain the true causal variant, whereby each variant is considered causal in turn and the predictions are normalised by that variant's posterior probability of causality.

------------------------------------------------------------------------

Abstract
--------

The primary goal of Genome Wide Association Studies (GWAS) is to better understand the biology of disease. GWAS have been successful in identifying thousands of associations with common and complex diseases, but these associations refer to physical genomic regions rather than specific causal variants responsible for the disease. Consequently, known GWAS association signals are used in follow-up studies, which continue to untangle the relationship between genetic variation and disease. Accurately localising the specific “causal variants” driving the association signals identified by GWAS is difficult due to correlations between SNPs. Standard practice in fine-mapping experiments is therefore to report a credible set of variants which is believed to contain the causal variant with some “coverage probability”.

We evaluated coverage probabilities of credible sets obtained using the dominant method in the field and found that the claimed coverage probabilities were systematically biased. In low power studies the true coverage was below that claimed, suggesting that researchers should add more variants to the credible set in order to achieve the required coverage. In high power studies, the true coverage was higher than that claimed, potentially allowing for higher resolution through the removal of variants from the set, whilst still attaining the required coverage. The algorithm to create credible sets contains an ordering step, which aims to make the set as small as possible, but which is not accounted for when estimating coverage probabilities. We showed that it is this ordering step that induces coverage bias and have developed a method to estimate this bias using rapid simulations based on the observed SNP correlation structure. We provide R software for our method which provides the user with an accurate coverage estimate of the causal variant in the credible set. Obtaining accurate coverage estimates will allow researchers to adjust their credible set to achieve the true desired coverage - for example narrowing down a 95% credible set of 10 variants to just 5 variants.

While technical, our result - that standard coverage estimates of the causal variant in GWAS are inaccurate and can be improved - impacts standard practice in genetic association studies. Improvement in the resolution of the fine-mapping experiment will enable more efficient expenditure of resources in the follow-up process of annotating the variants in the credible set to determine the implicated genes and pathways, helping to untangle the complex relationship between genetic variants and disease.

------------------------------------------------------------------------

[1] This refers to whether we include the null model of no genetic effect into the analyses. Note that in the standard Bayesian approach for fine-mapping, this is ignored such that the sum of the posterior probabilities for each variant is 1.
