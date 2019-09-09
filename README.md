
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
    (see [Corrected Coverage vignette](https://annahutch.github.io/corrcoverage/articles/corrected-coverage.html))
  - Find a new corrected credible set with the desired coverage of the
    true causal variant (see [Corrected Credible Set vignette](https://annahutch.github.io/corrcoverage/articles/New-Credible-Set.html))

We've strived to make our R package as easy to use as possible. 
Please see the flowchart below to decide which function is best to solve
your problem. The interactive version (click-to-functions) is [available
here](https://annahutch.github.io/PhD/package_flowchart.html).

![](https://annahutch.github.io/PhD/package_flowchart.svg)

-----

## Installation

We recommend that all users download the package straight from github, as this contains a more complete version of the package:

``` r
install.packages("devtools") # if not already installed
devtools::install_github("annahutch/corrcoverage")
```

Alternatively (and if using a solaris operating system) download straight from cran using:

```r
install.packages("corrcoverage")
```

-----

## Examples

For examples, please see the relevant vignettes.

1. The Corrected Coverage vignette
[here](https://annahutch.github.io/corrcoverage/articles/corrected-coverage.html)
should be read first. This shows readers how to use the `corrcoverage` R
package to get an accurate coverage estimate of the causal variant in a
credible set.

2. The Corrected Credible Set vignette
[here](https://annahutch.github.io/corrcoverage/articles/New-Credible-Set.html)
follows on from the Corrected Coverage vignette and shows readers how
the `corrcoverage` R package can be used to obtain a new credible set with
the desired coverage of the causal variant.

3. The Useful Info vignette
[here](https://annahutch.github.io/corrcoverage/articles/Useful-Info.html)
provides supplementary information about the usage of the package,
including information about other useful functions.

-----

In brief, the correction method involves simulating many credible sets
from the same system as the original and calculating what proportion of
these contain the true causal variant. Since the true causal variant is 
unknown, each variant is considered as causal in turn and the proportions 
are normalised by that variant’s posterior probability of causality.

-----

## Abstract

Genome Wide Association Studies (GWAS) have successfully identified thousands of loci associated with human diseases. Bayesian genetic fine-mapping studies aim to identify the specific causal variants within GWAS loci responsible for each association, reporting credible sets of plausible causal variants, which are interpreted as containing the causal variant with some "coverage probability".

Here, we investigate the coverage probabilities of credible sets through simulations and find that these are systematically biased. We present a method to re-estimate the coverage of credible sets using rapid simulations based on the observed, or estimated, SNP correlation structure, we call this the “corrected coverage estimate”. This is extended to find “corrected credible sets”, which are the smallest set of variants such that their corrected coverage estimate meets the target coverage. 

We use our method to improve the resolution of a fine-mapping study of type 1 diabetes. We found that in 27 out of 39 associated genomic regions our method could reduce the number of potentially causal variants to consider for follow-up, and found that none of the 95% or 99% credible sets required the inclusion of more variants – a pattern matched in simulations of well powered GWAS.

Crucially, our correction method requires only GWAS summary statistics and remains accurate when SNP correlations are estimated from a large reference panel. Using our method to improve the resolution of fine-mapping studies will enable more efficient expenditure of resources in the follow-up process of annotating the variants in the credible set to determine the implicated genes and pathways in human diseases. 

-----
