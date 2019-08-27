## Resubmission
This is a resubmission. In this version I have:

* Updated description field in DESCRIPTION file to include citations
  and remove incorrect quotation marks. 

* Updated authors of functions and correctly added these in the 
  Authors@R field with the appropriate roles.
  
* Replaced \dontrun{} with \donttest{} in long running (> 5 sec)
  examples.
  
* Removed any commented out code lines from examples.

* Removed code from vignettes that changes user's options
  (options(width = 3000)).

Test environments
-----------------

-   local OS X install, R 3.6.0
-   Ubuntu (on travis-ci, oldrel), R 3.5.3
-   Ubuntu (on travis-ci, release), R 3.6.1
-   Ubuntu (on travis-ci, devel), “R Under development (unstable)
    (2019-08-19 r77038)”
-   windows (using Win-Builder)

R CMD check results
-------------------

There were no ERRORs or WARNINGs or NOTEs.
