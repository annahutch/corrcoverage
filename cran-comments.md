## Resubmission
This is a resubmission. In this version I have:

* Updated title of package

* Removed simGWAS dependency

* Removed C++ zj_pp_alma function which was throwing errors on solaris check

## Test environments

-   local OS X install, R 3.6.0
-   Ubuntu (on travis-ci, oldrel), R 3.5.3
-   Ubuntu (on travis-ci, release), R 3.6.1
-   Ubuntu (on travis-ci, devel), “R Under development (unstable)
    (2019-08-19 r77038)”
-   Windows Server 2008

## R CMD check results

There were no ERRORs or WARNINGs or NOTEs.

## rhub check results

One NOTE:

checking installed package size ... NOTE
  installed size is  5.1Mb
  sub-directories of 1Mb or more:
      extdata   3.2Mb

This data is required for the vignettes. In the previous version submitted to cran, the simGWAS package was used in the vignettes to simulate GWAS summary stats. simGWAS has errors on cran so this data must be stored in the package directly.
