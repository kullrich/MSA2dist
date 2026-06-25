# MSA2dist <a href="https://www.bioconductor.org/packages/release/bioc/html/MSA2dist.html"><img src="man/figures/logo.png" align="right" height="160" /></a>

`MSA2dist` calculates pairwise distances between all sequences of a `DNAStringSet` or a `AAStringSet` using a custom score matrix and conducts codon based analysis. It uses scoring matrices to be used in these pairwise distance calculations which can be adapted to any scoring for DNA or AA characters. E.g. by using literal distances `MSA2dist` calculates pairwise `IUPAC` distances. `DNAStringSet` alignments can be analysed as codon alignments to look for synonymous and nonsynonymous substitutions (dN/dS) in a parallelised fashion using a variety of substitution models. Non-aligned coding sequences can be directly used to construct pairwise codon alignments (global/local) and calculate dN/dS without any external dependencies. In addition, `MSA2dist` provides population genetic analyses, including the calculation of nucleotide divergence between populations (Dxy) and genetic differentiation statistics (FST) from aligned sequence data.

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

## Supported dN/dS Models

Models used and implemented according to Li (1993) (via [seqinr](https://github.com/lbbe-software/seqinr)) and Nei & Gojobori (1986) (native implementation). In addition, the complete set of dN/dS estimation methods available in [KaKs_Calculator2](https://github.com/kullrich/kakscalculator2) has been ported and reimplemented in `MSA2dist` using [Rcpp](https://github.com/kullrich/MSA2dist/tree/devel/src), enabling efficient and dependency-free calculation of dN, dS, and dN/dS statistics directly within R.

| Model | Description |
|---------|-------------|
| `Li` | Li (1993) method |
| `NG86` | Nei & Gojobori (1986) method |
| `NG` | Nei & Gojobori method |
| `LWL` | Li-Wu-Luo method |
| `LPB` | Li-Pamilo-Bianchi method |
| `MLWL` | Modified Li-Wu-Luo method |
| `MLPB` | Modified Li-Pamilo-Bianchi method |
| `GY` | Goldman-Yang maximum-likelihood model |
| `YN` | Yang-Nielsen method |
| `MYN` | Modified Yang-Nielsen method |
| `MS` | Model Selection method |
| `MA` | Model Averaging method |
| `GNG` | Gamma-series Nei-Gojobori method |
| `GLWL` | Gamma-series Li-Wu-Luo method |
| `GLPB` | Gamma-series Li-Pamilo-Bianchi method |
| `GMLWL` | Gamma-series Modified Li-Wu-Luo method |
| `GMLPB` | Gamma-series Modified Li-Pamilo-Bianchi method |
| `GYN` | Gamma-series Yang-Nielsen method |
| `GMYN` | Gamma-series Modified Yang-Nielsen method |

These models differ in their assumptions regarding codon frequencies, transition/transversion bias, unequal substitution rates among sites, and rate heterogeneity. This allows users to select simple counting-based approaches for rapid screening or more sophisticated maximum-likelihood and gamma-corrected methods for evolutionary analyses.

## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://github.com/kullrich/MSA2dist/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.

