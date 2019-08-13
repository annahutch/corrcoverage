
<!-- README.md is generated from README.Rmd. Please edit that file -->

# corrcoverage <img src="man/figures/logo.png" align="right" />

[![Coverage
status](https://codecov.io/gh/annahutch/corrcoverage/branch/master/graph/badge.svg)](https://codecov.io/github/annahutch/corrcoverage?branch=master)
[![Build
Status](https://travis-ci.org/annahutch/corrcoverage.svg?branch=master)](https://travis-ci.org/annahutch/corrcoverage)

Webpage: <https://annahutch.github.io/corrcoverage/>

The `corrcoverage` R package uses a computationally efficient algorithm
to find accurate coverage estimates of the causal variant in credible
sets of prioritised variants from single causal variant genetic fine mapping ([Maller et
al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/23104008),
[Wakefield, 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359).).

The package only requires GWAS summary statistics and can be used to:

  - Perform single causal variant Bayesian fine-mapping
  - Estimate the true genetic effect at the causal variant (see `est_mu`
    function)
  - Obtain an accurate coverage estimate that the causal variant 
    is contained within a credible set, the "corrected coverage estimate"
    (see `corrcov` function and "Corrected Coverage"" vignette)
  - Find a new "corrected" credible set with the desired coverage of the
    true causal variant (see `corrected_cs` function and 
    "Corrected Credible Set"" vignette)

We've strived to make our R package as easy to use as possible. 
Please see the flowchart below to decide which function is best to solve
your problem. The interactive version (click-to-functions) is [available
here](https://annahutch.github.io/PhD/package_flowchart.html).

![](https://annahutch.github.io/PhD/package_flowchart.svg)

-----

## Installation

You can install the released version of `corrcoverage` from
[github](https://github.com/) with:

``` r
install.packages("devtools") # if not already installed
devtools::install_github("annahutch/corrcoverage")
```

-----

## Examples

For examples, please see the relevant vignettes.

1. The "Corrected Coverage" vignette
[here](https://annahutch.github.io/corrcoverage/articles/Corrected-Coverage.html)
should be read first. This shows readers how to use the `corrcoverage` R
package to get an accurate coverage estimate of the causal variant in a
credible set.

2. The "Corrected Credible Set" vignette
[here](https://annahutch.github.io/corrcoverage/articles/New-Credible-Set.html)
follows on from the "Corrected Coverage"" vignette and shows readers how
the `corrcoverage` R package can be used to obtain a new credible set with
the desired coverage of the causal variant.

3. The ‘Useful Info’ vignette
[here](https://annahutch.github.io/corrcoverage/articles/Useful-Info.html)
provides supplementary information about the usage of the package,
including information about other useful functions.

-----

In brief, the correction method involves simulating many credible sets
from "the same system as the original" and calculating what proportion of
these contain the true causal variant. Since the true causal variant is 
unknown, each variant is considered as causal in turn and the proportions 
are normalised by that variant’s posterior probability of causality.

-----

## Abstract

Genome Wide Association Studies (GWAS) have successfully identified thousands of loci associated with human diseases. Bayesian genetic fine-mapping studies aim to identify the specific causal variants within GWAS loci responsible for the disease, reporting credible sets of plausible causal variants, which are interpreted as containing the causal variant with some ‘coverage probability’.

Here, we investigate the coverage probabilities of credible sets through simulations and find them to be systematically biased. We present a method to re-estimate the coverage of credible sets using rapid simulations based on the observed, or estimated, SNP correlation structure, we call this re-estimation the “corrected coverage estimate”. This is extended to find “corrected credible sets”, which are the smallest set of variants such that the corrected coverage probability is accurate. 

We show that our method can be used to improve the resolution of fine-mapping studies using empirical type 1 diabetes data. Specifically, we found that in 29 out of 39 associated genomic regions, our method could be used to reduce the number of potentially causal variants to consider for follow-up – leading to more powerful and reliable follow up studies. 

Crucially, our correction method requires only genetic association summary test statistics and remains accurate when SNP correlations are estimated from a large reference panel. Using our method to improve the resolution of fine-mapping experiments will enable more efficient expenditure of resources in the follow-up process of annotating the variants in the credible set to determine the implicated genes and pathways. 

-----
