# MSA2dist 1.7.3 (2024-03-14)

BUG FIXES

* fixed dnastring2kaks to return correct orientation
* fixed pal2nal to use gap_pos-n_i_codons_added-1

# MSA2dist 1.7.2 (2024-03-14)

BUG FIXES

* fixed pal2nal to cover individual gaps longer than 1

# MSA2dist 1.7.1 (2024-03-13)

NEW FEATURES

* added pal2nal
* added cdsstring2codonaln
* added indices2kaks
* added KaKs Calculator 2.0 example to vignette

CHANGES

* changed dnastring2kaks to use cdsstring2codonaln function

# MSA2dist 1.5.2 (2023-05-24)

NEW FEATURES

* added aa2selfscore

# MSA2dist 1.5.1 (2023-05-22)

NEW FEATURES

* additional option to use asymmetric score matrix (symmetric=FALSE)

CHANGES
    
* changed rcpp_distSTRING.cpp
* changed rcpp_pairwiseDeletionAA.cpp
* changed rcpp_pairwiseDeletionDNA.cpp
* changed aastring2dist.R
* changed dnastring2dist.R

# MSA2dist 1.3.1 (2022-11-08)

NEW FEATURES

* additional Genetic Codes into base.cpp

# MSA2dist 1.1.7 (2022-10-05)

CHANGES

* changed licence into GPL-3 to account for KaKs Calculator 2.0 licence

# MSA2dist 1.1.6 (2022-09-22)

CHANGES

* changed test-rcpp_KaKs

# MSA2dist 1.1.5 (2022-09-19)

CHANGES

* changed rcpp_KaKs data access from data.frame to list

# MSA2dist 1.1.4 (2022-07-16)

BUG FIXES

* fixed rcpp Warnings

# MSA2dist 1.1.3 (2022-07-11)

BUG FIXES

* fixed rcpp indentation Warnings

# MSA2dist 1.1.2 (2022-07-05)

NEW FEATURES

* Added KaKs Calculator 2.0 models as rcpp implementation

BUG FIXES

* fixed some typos

# MSA2dist 1.1.1 (2022-06-30)

NEW FEATURES

* Added cds2codonaln

# MSA2dist 0.99.3 (2022-03-18)

BUG FIXES

* Fix RcppThread::LdFlags warning

# MSA2dist 0.99.2 (2022-01-28)

NEW FEATURES

* Added RcppThread::ProgressBar

# MSA2dist 0.99.1 (2022-01-27)

CHANGES

* Changed version number into 0.99.2
* Changed URL links in DESCRIPTION

# MSA2dist 0.99.0 (2021-12-22)

NEW FEATURES

* Submitted to Bioconductor

CHANGES

* Changed version number into 0.99.1
* Changed name from distSTRING into MSA2dist

