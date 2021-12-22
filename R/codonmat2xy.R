#' @title codonmat2xy
#' @name codonmat2xy
#' @description This function calculates average behavior of each codon for all
#' pairwise comparisons for indels, syn, and nonsyn mutations according to
#' \emph{Nei and Gojobori (1986)}.
#' @param codonmat \code{codon matrix} obtained via
#' \code{\link[MSA2dist]{dnastring2codonmat}} [mandatory]
#' @param threads number of parallel threads [default: 1]
#' @return A \code{data.frame} object with the following components:\cr
#' \code{Codon} Codon index\cr
#' \code{n} number of comparison\cr
#' \code{SynSum} Sum of syn\cr
#' \code{NonSynSum} Sum of nonsyn\cr
#' \code{IndelSum} Sum of indels\cr
#' \code{SynMean} average syn per codon\cr
#' \code{NonSynMean} average nonsyn per codon\cr
#' \code{IndelMean} average indels per codon\cr
#' \code{CumSumSynMean} cumulative average syn per codon\cr
#' \code{CumSumNonSynMean} cumulative average nonsyn per codon\cr
#' \code{CumSumIndelMean} cumulative indels per codon\cr
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom seqinr kaks
#' @importFrom tidyr %>% unite
#' @importFrom tibble add_column
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom parallel makeForkCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr group_by filter count left_join summarise mutate
#' @importFrom rlang .data
#' @seealso \code{\link[MSA2dist]{dnastring2codonmat}}
#' \code{\link[MSA2dist]{codonmat2pnps}}
#' \code{\link[MSA2dist]{dnastring2kaks}}
#' \code{\link[seqinr]{kaks}}
#' @references Nei and Gojobori. (1986) Simple methods for estimating the
#' numbers of synonymous and nonsynonymous nucleotide substitutions.
#' \emph{Mol. Biol. Evol.}, \bold{3(5)}, 418-426.
#' @references Ganeshan et al. (1997) Human immunodeficiency virus type 1
#' genetic evolution in children with different rates of development of disease.
#' \emph{J. Virology.} \bold{71(1)}, 663-677.
#' @references Yang et al. (2000) Codon-substitution models for heterogeneous
#' selection pressure at amino acid sites. \emph{Genetics.}
#' \bold{155(1)}, 431-449.
#' @examples
#' ## load example sequence data
#' data("hiv", package="MSA2dist")
#' #codonmat2xy(dnastring2codonmat(hiv))
#' hiv |> dnastring2codonmat() |> codonmat2xy()
#' #codonmat2xy(dnastring2codonmat(hiv), threads=2)
#' hiv |> dnastring2codonmat() |> codonmat2xy(threads=2)
#' @export codonmat2xy
#' @author Kristian K Ullrich

codonmat2xy <- function(codonmat, threads = 1){
    if(.Platform$OS.type == "windows"){
        cl <- parallel::makeCluster(threads)
    }
    if(.Platform$OS.type != "windows"){
        cl <- parallel::makeForkCluster(threads)
    }
    doParallel::registerDoParallel(cl)
    i <- NULL
    j <- NULL
    k <- NULL
    OUT <- foreach::foreach(i = seq(from = 1, to = nrow(codonmat)),
        .combine=rbind, .packages = c('foreach')) %dopar% {
        codon_i <- table(codonmat[i,])
####
#        foreach::foreach(j = seq(from = 1, to = ncol(codonmat) - 1),
#            .combine=rbind) %do% {
#            foreach::foreach(k = seq(from = j + 1, to = ncol(codonmat)),
#                .combine=rbind) %do% {
#        c(setNames(i, "Codon"),
#            setNames(j, "Comp1"),
#            setNames(k, "Comp2"),
#            setNames(MSA2dist::compareCodons(codonmat[i, j],
#            codonmat[i, k]), c("syn", "nonsyn", "indel")))
####
        foreach::foreach(j = seq(from = 1, to = length(codon_i) ),
            .combine=rbind) %do% {
            foreach::foreach(k = seq(from = j, to = length(codon_i)),
                .combine=rbind) %do% {
                if(j == k){
                    d <- do.call("rbind", rep(list(c(setNames(i, "Codon"),
                        setNames(j, "Comp1"),
                        setNames(k, "Comp2"),
                        setNames(MSA2dist::compareCodons(
                            names(codon_i)[j],
                            names(codon_i)[k]),
                            c("syn", "nonsyn", "indel")))),
                    (codon_i[j] * (codon_i[k] - 1)) / 2))
                }
                if(j != k){
                    d <- do.call("rbind", rep(list(c(setNames(i, "Codon"),
                        setNames(j, "Comp1"),
                        setNames(k, "Comp2"),
                        setNames(MSA2dist::compareCodons(
                            names(codon_i)[j],
                            names(codon_i)[k]),
                            c("syn", "nonsyn", "indel")))),
                    codon_i[j] * codon_i[k]))
                }
                d
            }
        }
    }
    parallel::stopCluster(cl)
    OUT <- as.data.frame(OUT)
    OUT.NAs <- OUT %>% dplyr::group_by(.data$Codon) %>%
        dplyr::filter(!is.na(.data$syn)) %>%
        dplyr::count(.data$Codon)
    OUT.SynSum <- OUT %>% dplyr::group_by(.data$Codon) %>%
        dplyr::summarise(SynSum = sum(.data$syn, na.rm = TRUE))
    OUT.NonSynSum <- OUT %>% dplyr::group_by(.data$Codon) %>%
        dplyr::summarise(NonSynSum = sum(.data$nonsyn, na.rm = TRUE))
    OUT.IndelSum <- OUT %>% dplyr::group_by(.data$Codon) %>%
        dplyr::summarise(IndelSum = sum(.data$indel, na.rm = TRUE))
    OUT.join <- dplyr::left_join(OUT.NAs, OUT.SynSum) %>%
        dplyr::left_join(OUT.NonSynSum) %>% dplyr::left_join(OUT.IndelSum)
    OUT.xy <- OUT.join %>% dplyr::mutate(SynMean = .data$SynSum/.data$n,
        NonSynMean = .data$NonSynSum/.data$n,
        IndelMean = .data$IndelSum/.data$n)
    OUT.xy <- OUT.xy %>%
        tibble::add_column(CumSumSynMean = cumsum(OUT.xy$SynMean),
        CumSumNonSynMean = cumsum(OUT.xy$NonSynMean),
        CumSumIndelMean = cumsum(OUT.xy$IndelMean))
    return(OUT.xy)
}
