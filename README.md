# MSA2dist <a href="https://www.bioconductor.org/packages/release/bioc/html/MSA2dist.html"><img src="man/figures/logo.png" align="right" height="160" /></a>

`MSA2dist` calculates pairwise distances between all sequences of a `DNAStringSet` or a `AAStringSet` using a custom score matrix and conducts codon based analysis.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `MSA2dist` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("MSA2dist")
```

And the development version from
[GitHub](https://github.com/kullrich/MSA2dist) with:

``` r
BiocManager::install("kullrich/MSA2dist")
```

## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://github.com/kullrich/MSA2dist/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.

