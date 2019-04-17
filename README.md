
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->
[![Travis build
status](https://travis-ci.org/tidyverse/dplyr.svg?branch=master)](https://travis-ci.org/annahutch/corrcoverage)
<!-- badges: end -->

corrcoverage
============

Webpage: <https://annahutch.github.io/corrcoverage/>

The `corrcoverage` R package uses a computationally efficient algorithm to find accurate coverage estimates of the causal variant in credible sets obtained using the Bayesian approach for fine-mapping ([Maller et al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/23104008), [Wakefield, 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359).).

The package only requires GWAS summary statistics and can be used to:

-   Perform Bayesian fine-mapping
-   Estimate the true genetic effect at the causal variant
-   Obtain an accurate coverage estimate of the causal variant in a credible set (corrected coverage estimate)
-   Find the smallest set of variants such that the coverage estimate exceeds some user-defined threshold (corrected credible set)

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

For examples, please see the vignettes on the webpage [here](http://annahutch.github.io/corrcoverage/articles/my-vignette.html).

The 'Corrected Coverage' vignette [here](https://annahutch.github.io/corrcoverage/articles/my-vignette.html) should be read first. This shows readers how to use the `corrcoverage` package to get an accurate coverage estimate of the causal variant in a credible set.

The 'New Credible Set' vignette [here](https://annahutch.github.io/corrcoverage/articles/New-Credible-Set.html) follows on from the 'Corrected Coverage' vignette and shows readers how the `corrcoverage` package can be used to obtain a new credible set with the desired coverage of the causal variant.

------------------------------------------------------------------------

Conversion Functions
--------------------

The Bayesian method for fine-mapping involves finding the posterior probability of causality for each SNP, before sorting these into descending order and adding variants to a 'credible set' until the combined posterior probabilities of these SNPs exceed some threshold. The supplementary text of Maller's paper (available [here](https://media.nature.com/original/nature-assets/ng/journal/v44/n12/extref/ng.2435-S1.pdf)) shows that these posterior probabilities are normalised Bayes factors.

Asymptotic Bayes factors ([Wakefield, 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359)) are commonly used in genetic association studies as these only require the specification of *Z*-scores (or equivalently the effect size coefficients, *β*, and their standard errors, *V*), the standard errors of the effect sizes (*V*) and the prior variance of the estimated effect size (*W*<sup>2</sup>), thus only requiring summary data from genetic association studies.

Consequently, the `corrcoverage` package contains functions for converting between *P*-values, *Z*-scores, asymptotic Bayes factors (ABFs) and posterior probabilities of causality (PPs). The following table shows what input these conversion functions require and what output they produce. The 'include null model' column is for whether the null model of no genetic effect is included in the calculation (PPs obtained using the standard Bayesian approach ignore this).

| Function      | Include null model? | Input                     | Output                  |
|---------------|---------------------|---------------------------|-------------------------|
| `approx.bf.p` | YES                 | *p*-values                | log(ABF)                |
| `pvals_pp`    | YES                 | *p*-values                | Posterior Probabilities |
| `z0_pp`       | YES                 | Marginal *Z*-scores       | Posterior Probabilities |
| `ppfunc`      | NO                  | Marginal *Z*-scores       | Posterior Probabilities |
| `z0func.mat`  | NO                  | Marginal *Z*-score matrix | Posterior Probabilities |

Functions are also provided to simulate marginal *Z*-scores from joint *Z*-scores. The joint *Z*-scores are all 0, except at the causal variant where it is the 'true effect', *μ*.

*μ* can be estimated using the `est_mu` function which required sample sizes, *Z*-scores and minor allele frequencies.

`z_sim` simulates marginal *Z*-scores from joint *Z*-scores, whilst `zj_pp` can be used to simulate posterior probability systems from a joint *Z*-score vector. These functions first calculate *E*(*Z*<sub>*m*</sub>),
*E*(*Z*<sub>*m*</sub>)=*Z*<sub>*j*</sub> × *Σ*
 where *Σ* is the correlation matrix between SNPs.

We can then simulate more *Z*-score systems from a multivariate normal distribution with mean *E*(*Z*<sub>*m*</sub>) and variance *Σ*. This is a key step in our corrected coverage method.

------------------------------------------------------------------------

In brief, the correction method includes simulating many credible sets from 'the same system as the original' and finding what proportion of these contain the true causal variant, whereby each variant is considered causal in turn and the predictions are normalised by that variant's posterior probability of causality.

------------------------------------------------------------------------

Abstract
--------

The primary goal of Genome Wide Association Studies (GWAS) is to better understand the biology of disease. GWAS have been successful in identifying thousands of associations with common and complex diseases, but these associations refer to physical genomic regions rather than specific causal variants responsible for the disease. Consequently, known GWAS association signals are used in follow-up studies, which continue to untangle the relationship between genetic variation and disease. Accurately localising the specific “causal variants” driving the association signals identified by GWAS is difficult due to correlations between SNPs. Standard practice in fine-mapping experiments is therefore to report a credible set of variants which is believed to contain the causal variant with some “coverage probability”.

We evaluated coverage probabilities of credible sets obtained using the dominant method in the field and found that the claimed coverage probabilities were systematically biased. In low power studies the true coverage was below that claimed, suggesting that researchers should add more variants to the credible set in order to achieve the required coverage. In high power studies, the true coverage was higher than that claimed, potentially allowing for higher resolution through the removal of variants from the set, whilst still attaining the required coverage. The algorithm to create credible sets contains an ordering step, which aims to make the set as small as possible, but which is not accounted for when estimating coverage probabilities. We showed that it is this ordering step that induces coverage bias and have developed a method to estimate this bias using rapid simulations based on the observed SNP correlation structure. We provide R software for our method which provides the user with an accurate coverage estimate of the causal variant in the credible set. Obtaining accurate coverage estimates will allow researchers to adjust their credible set to achieve the true desired coverage - for example narrowing down a 95% credible set of 10 variants to just 5 variants.

While technical, our result - that standard coverage estimates of the causal variant in GWAS are inaccurate and can be improved - impacts standard practice in genetic association studies. Improvement in the resolution of the fine-mapping experiment will enable more efficient expenditure of resources in the follow-up process of annotating the variants in the credible set to determine the implicated genes and pathways, helping to untangle the complex relationship between genetic variants and disease.

------------------------------------------------------------------------
