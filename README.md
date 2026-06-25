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

## Examples

### Computing pairwise dN/dS using parallel threads

```
## load example sequence data
data("hiv", package="MSA2dist")
hiv.dNdS <- hiv |> dnastring2kaks(model="MYN", threads=2)
```

```
> tibble(hiv.dNdS)
# A tibble: 78 × 25
   Comp1 Comp2 seq1   seq2   Method Ka    Ks    `Ka/Ks` `P-Value(Fisher)` Length
   <chr> <chr> <chr>  <chr>  <chr>  <chr> <chr> <chr>   <chr>             <chr> 
 1 1     2     U68496 U68497 MYN    0.02… 0.09… 0.2672… 0.104824          273   
 2 1     3     U68496 U68498 MYN    0.08… 0.04… 1.85833 0.665602          273   
 3 1     4     U68496 U68499 MYN    0.08… 0.05… 1.66318 0.645824          273   
 4 1     5     U68496 U68500 MYN    0.11… 0.11… 1.00159 0.96729           273   
 5 1     6     U68496 U68501 MYN    0.10… 0.07… 1.35518 0.731834          273   
 6 1     7     U68496 U68502 MYN    0.14… 0.17… 0.8292… 0.600487          273   
 7 1     8     U68496 U68503 MYN    0.13… 0.09… 1.41146 0.589764          273   
 8 1     9     U68496 U68504 MYN    0.11… 0.23… 0.4656… 0.0955289         273   
 9 1     10    U68496 U68505 MYN    0.10… 0.27… 0.4019… 0.116874          273   
10 1     11    U68496 U68506 MYN    0.11… 0.11… 1.05831 0.944884          273   
# ℹ 68 more rows
# ℹ 15 more variables: `S-Sites` <chr>, `N-Sites` <chr>,
#   `Fold-Sites(0:2:4)` <chr>, Substitutions <chr>, `S-Substitutions` <chr>,
#   `N-Substitutions` <chr>, `Fold-S-Substitutions(0:2:4)` <chr>,
#   `Fold-N-Substitutions(0:2:4)` <chr>, `Divergence-Time` <chr>,
#   `Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)` <chr>,
#   `GC(1:2:3)` <chr>, `ML-Score` <chr>, AICc <chr>, `Akaike-Weight` <chr>, …
# ℹ Use `print(n = ...)` to see more rows
```

### Computing pairwise dN/dS with automatic pairwise codon alignments from unaligned sequences (option `isMSA = FALSE`)

```
## define three unaligned cds sequences
cds1 <- Biostrings::DNAString("ATGCAACATTGC")
cds2 <- Biostrings::DNAString("ATGCATTGC")
cds3 <- Biostrings::DNAString("ATGCAATGC")
cds_sequences <- Biostrings::DNAStringSet(list(cds1, cds2, cds3))
names(cds_sequences) <- c("cds1", "cds2", "cds3")
cds_sequences |> dnastring2kaks(model="Li", isMSA=FALSE)
```

Alignment parameter can be adapted (see `?cds2codonaln` for options)

```
type = "global",
substitutionMatrix = "BLOSUM62",
gapOpening = 10,
gapExtension = 0.5,
remove.gaps = FALSE,
```

### Computing dxy and fst for all population combinations of a ‘DNAStringSet’

```
## load example sequence data
data("iupac", package="MSA2dist")
poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
    GER = grep("Mmd.GER", names(iupac)),
    IRA = grep("Mmd.IRA", names(iupac)),
    AFG = grep("Mmm.AFG", names(iupac)))
iupac <- iupac |> addpop2string(poplist)
dxy <- iupac |> dnastring2dxy(model="IUPAC", pop=poplist)
```

```
> dxy$dxy
            FRA         GER         IRA          AFG
FRA 0.002304687 0.004015625 0.005171875 0.0122536740
GER 0.004015625 0.003835937 0.005703125 0.0123789245
IRA 0.005171875 0.005703125 0.006320313 0.0118779225
AFG 0.012253674 0.012378925 0.011877923 0.0008628368

> dxy$fst
          FRA       GER       IRA       AFG
FRA 0.0000000 0.2354086 0.1661631 0.8707521
GER 0.2354086 0.0000000 0.1095890 0.8102107
IRA 0.1661631 0.1095890 0.0000000 0.6976260
AFG 0.8707521 0.8102107 0.6976260 0.0000000
```

## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://github.com/kullrich/MSA2dist/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.

