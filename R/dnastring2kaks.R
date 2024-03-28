#' @title dnastring2kaks
#' @name dnastring2kaks
#' @description This function calculates Ka/Ks (pN/pS)
#' for all combinations of a \code{DNAStringSet}.
#' If the sequences in the \code{DNAStringSet} are not a multiple-sequence
#' alignment, pairwise codon alignments can be calculated on the fly.
#' Models used and implemented according to
#' \emph{Li (1993)} (using \code{seqinr}) or
#' \emph{Nei and Gojobori (1986)} (own implementation) or models from
#' \code{KaKs_Calculator2} ported to \code{MSA2dist} with \code{Rcpp}.
#' @param cds \code{DNAStringSet} coding sequence alignment [mandatory]
#' @param model specify codon model either "Li" or "NG86" or
#' one of \code{KaKs_Calculator2} model "NG", "LWL", "LPB", "MLWL", "MLPB",
#' "GY", "YN", "MYN", "MS", "MA", "GNG", "GLWL", "GLPB", "GMLWL", "GMLPB",
#' "GYN", "GMYN" [default: Li]
#' @param threads number of parallel threads [default: 1]
#' @param isMSA cds \code{DNAStringSet} represents MSA [default: TRUE]
#' @param sgc standard genetic code (for KaKs Calculator models)
#' [default: 1]
#' @param verbose verbosity (for KaKs Calculator models) [default: FALSE]
#' @param ... other codon alignment parameters
#' @return A \code{data.frame} of \code{KaKs} values
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom seqinr kaks
#' @importFrom parallel makeForkCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom tidyr %>% as_tibble pivot_longer
#' @importFrom dplyr slice left_join
#' @importFrom tibble column_to_rownames add_column
#' @importFrom stringr str_split
#' @seealso \code{\link[seqinr]{kaks}}
#' @references "MS/MA/GNG/GLWL/GLPB/GMLWL/GMLPB/GYN:" Wang et al. (2010)
#' KaKs_Calculator 2.0: a toolkit incorporating
#' gamma-series methods and sliding window strategies.\emph{Genomics,
#' proteomics & bioinformatics.} \bold{8(1)}, 77-80.
#' @references "Li/LWL:" Li et al. (1985) A new method for estimating synonymous
#' and nonsynonymous rates of nucleotide substitution considering the relative
#' likelihood of nucleotide and codon changes. \emph{Mol. Biol. Evol.},
#' \bold{2(2)}, 150-174.
#' @references "Li/LPB:" Li (1993). Unbiased estimation of the rates of
#' synonymous and nonsynonymous substitution. Journal of molecular evolution,
#' 36(1), pp.96-99.
#' @references "NG86/NG:" Nei and Gojobori. (1986) Simple methods for estimating
#' the numbers of synonymous and nonsynonymous nucleotide substitutions.
#' \emph{Mol. Biol. Evol.}, \bold{3(5)}, 418-426.
#' @references "LPB:" Pamilo and Bianchi. (1993) Evolution of the Zfx and Zfy
#' genes: Rates and interdependence between genes. \emph{Mol. Biol. Evol.},
#' \bold{10}, 271-281.
#' @references "MLWL/MLPB:" Tzeng et al. (2004). Comparison of three methods for
#' estimating rates of synonymous and nonsynonymous nucleotide substitutions.
#' \emph{Mol. Biol. Evol.}, \bold{21(12)}, 2290-2298.
#' @references "GY:" Goldman and Yang (1994). A codon-based model of nucleotide
#' substitution for protein-coding DNA sequences. \emph{Mol. Biol. Evol.},
#' \bold{11(5)} 725-736.
#' @references "YN:" Yang et al. (2000) Codon-substitution models for
#' heterogeneous selection pressure at amino acid sites. \emph{Genetics.}
#' \bold{155(1)}, 431-449.
#' @references "MYN:" Zhang et al. (2006). Computing Ka and Ks with a
#' consideration of unequal transitional substitutions.
#' \emph{BMC evolutionary biology}, \bold{6(1)}, 1-10.
#' @references "data(hiv):" Ganeshan et al. (1997) Human immunodeficiency virus
#' type 1 genetic evolution in children with different rates of development of
#' disease. \emph{J. Virology.} \bold{71(1)}, 663-677.
#' @references Wang et al. (2009). gamma-MYN: a new algorithm for estimating Ka
#' and Ks with consideration of variable substitution rates.
#' \emph{Biology Direct}, \bold{4(1)}, 1-18.
#' @examples
#' ## load example sequence data
#' data("hiv", package="MSA2dist")
#' #dnastring2kaks(hiv, model="Li")
#' hiv |> dnastring2kaks(model="Li")
#' #dnastring2kaks(hiv, model="NG86")
#' hiv |> dnastring2kaks(model="NG86")
#' #dnastring2kaks(hiv, model="NG86", threads=2)
#' hiv |> dnastring2kaks(model="NG86", threads=2)
#'
#' ## define three unaligned cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATGCATTGC")
#' cds3 <- Biostrings::DNAString("ATGCAATGC")
#' cds_sequences <- Biostrings::DNAStringSet(list(cds1, cds2, cds3))
#' names(cds_sequences) <- c("cds1", "cds2", "cds3")
#' ## set isMSA to FALSE to automatically create pairwise codon alignments
#' #dnastring2kaks(cds_sequences, model="Li", isMSA=FALSE)
#' cds_sequences |> dnastring2kaks(model="Li", isMSA=FALSE)
#' @export dnastring2kaks
#' @author Kristian K Ullrich

dnastring2kaks <- function(cds,
    model = "Li",
    threads = 1,
    isMSA = TRUE,
    sgc = "1",
    verbose = FALSE,
    ...){
    stopifnot("Error: input needs to be a DNAStringSet"=
        methods::is(cds, "DNAStringSet"))
    local_cds2aa <- function(...,
        type,
        substitutionMatrix,
        gapOpening,
        gapExtension,
        remove.gaps) {MSA2dist::cds2aa(...)}
    local_cdsstring2codonaln <- function(...,
        shorten,
        frame,
        framelist,
        genetic.code,
        return.cds) {MSA2dist::cdsstring2codonaln(...)}
    Comp1 <- FALSE
    Comp2 <- FALSE
    seq1 <- FALSE
    seq2 <- FALSE
    cds.names <- names(cds)
    cds.names <- gsub(" ", "_",
        gsub(" $", "", gsub("\\s+", " ", cds.names)))
    names(cds) <- cds.names
    stopifnot("Error: either choose model 'Li' or 'NG86'
        or KaKs_Calculator2 model"=
        model %in% c("Li", "NG86",
            "NG", "LWL", "LPB", "MLWL", "MLPB", "GY", "YN",
            "MYN", "MS", "MA", "GNG", "GLWL", "GLPB", "GMLWL",
            "GMLPB", "GYN", "GMYN"))
    if(model == "Li"){
        if("Comp1" %in% names(cds)){
            names(cds)[which(names(cds) == "Comp1")] <- "_Comp1"
            Comp1 <- TRUE
        }
        if("Comp2" %in% names(cds)){
            names(cds)[which(names(cds) == "Comp2")] <- "_Comp2"
            Comp2 <- TRUE
        }
        if("seq1" %in% names(cds)){
            names(cds)[which(names(cds) == "seq1")] <- "_seq1"
            seq1 <- TRUE
        }
        if("seq2" %in% names(cds)){
            names(cds)[which(names(cds) == "seq2")] <- "_seq2"
            seq2 <- TRUE
        }
        if(isMSA){
            OUT <- seqinr::kaks(dnastring2aln(cds))
            OUT.ka <- as.matrix(OUT$ka)
            OUT.ks <- as.matrix(OUT$ks)
            OUT.vka <- as.matrix(OUT$vka)
            OUT.vks <- as.matrix(OUT$vks)
            OUT.ka <- OUT.ka %>% tidyr::as_tibble(rownames = "seq1") %>%
                tidyr::pivot_longer(colnames(OUT.ka),
                names_to = "seq2", values_to = "ka") %>%
                dplyr::slice(uptriidx(ncol(OUT.ka)))
            OUT.ks <- OUT.ks %>% tidyr::as_tibble(rownames = "seq1") %>%
                tidyr::pivot_longer(colnames(OUT.ks),
                names_to = "seq2", values_to = "ks") %>%
                dplyr::slice(uptriidx(ncol(OUT.ks)))
            OUT.vka <- OUT.vka %>% tidyr::as_tibble(rownames = "seq1") %>%
                tidyr::pivot_longer(colnames(OUT.vka),
                names_to = "seq2", values_to = "vka") %>%
                dplyr::slice(uptriidx(ncol(OUT.vka)))
            OUT.vks <- OUT.vks %>% tidyr::as_tibble(rownames = "seq1") %>%
                tidyr::pivot_longer(colnames(OUT.vks),
                names_to = "seq2", values_to = "vks") %>%
                dplyr::slice(uptriidx(ncol(OUT.vks)))
            OUT <- OUT.ka %>% dplyr::left_join(OUT.ks) %>%
                dplyr::left_join(OUT.vka) %>% dplyr::left_join(OUT.vks)
            OUT <- as.data.frame(cbind(Comp1=match(OUT$seq1, cds.names),
                Comp2=match(OUT$seq2, cds.names), OUT))
            if(Comp1){
                OUT$seq1 <- gsub("_Comp1", "Comp1", OUT$seq1)
                OUT$seq2 <- gsub("_Comp1", "Comp1", OUT$seq2)
            }
            if(Comp2){
                OUT$seq1 <- gsub("_Comp2", "Comp2", OUT$seq1)
                OUT$seq2 <- gsub("_Comp2", "Comp2", OUT$seq2)
            }
            if(seq1){
                OUT$seq1 <- gsub("_seq1", "seq1", OUT$seq1)
                OUT$seq2 <- gsub("_seq1", "seq1", OUT$seq2)
            }
            if(seq2){
                OUT$seq1 <- gsub("_seq2", "seq2", OUT$seq1)
                OUT$seq2 <- gsub("_seq2", "seq2", OUT$seq2)
            }
            OUT[["Ka"]] <- OUT[["ka"]]
            OUT[["Ks"]] <- OUT[["ks"]]
            OUT[["Ka/Ks"]] <- as.numeric(OUT[["ka"]])/
                as.numeric(OUT[["ks"]])
            attr(OUT, "model") <- "Li"
            attr(OUT, "align") <- "FALSE"
            attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT)
        } else {
            if(.Platform$OS.type == "windows"){
                cl <- parallel::makeCluster(threads)
            }
            if(.Platform$OS.type != "windows"){
                cl <- parallel::makeForkCluster(threads)
            }
            doParallel::registerDoParallel(cl)
            i <- NULL
            j <- NULL
            aa <- local_cds2aa(cds, ...)
            cds <- local_cds2aa(cds, return.cds=TRUE, ...)
            OUT <- foreach(i = seq(from = 1, to = length(cds) - 1),
                .combine=rbind, .packages = c('foreach')) %dopar% {
                foreach(j = seq(from = i + 1, to = length(cds)),
                    .combine=rbind) %do% {
                    c(setNames(i, "Comp1"),
                    setNames(j, "Comp2"),
                    setNames(cds.names[i], "seq1"),
                    setNames(cds.names[j], "seq2"),
                    unlist(seqinr::kaks(MSA2dist::dnastring2aln(
                        local_cdsstring2codonaln(cds[c(i, j)], aa[c(i, j)],
                            ...)))))
                }
            }
            parallel::stopCluster(cl)
            if(is.null(dim(OUT))){
                OUT <- as.data.frame(t(OUT))
            } else {
                OUT <- as.data.frame(OUT)
            }
            OUT[["Ka"]] <- OUT[["ka"]]
            OUT[["Ks"]] <- OUT[["ks"]]
            OUT[["Ka/Ks"]] <- as.numeric(OUT[["ka"]])/
                as.numeric(OUT[["ks"]])
            attr(OUT, "model") <- "Li"
            attr(OUT, "align") <- "TRUE"
            attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT)
        }
    }
    if(model == "NG86"){
        if(isMSA){
            if(.Platform$OS.type == "windows"){
                cl <- parallel::makeCluster(threads)
            }
            if(.Platform$OS.type != "windows"){
                cl <- parallel::makeForkCluster(threads)
            }
            doParallel::registerDoParallel(cl)
            i <- NULL
            j <- NULL
            codonmat <- MSA2dist::dnastring2codonmat(cds)
            OUT <- foreach(i = seq(from = 1, to = ncol(codonmat) - 1),
                .combine=rbind, .packages = c('foreach')) %dopar% {
                foreach(j = seq(from = i + 1, to = ncol(codonmat)),
                    .combine=rbind) %do% {
                    c(setNames(i, "Comp1"),
                    setNames(j, "Comp2"),
                    MSA2dist::codonmat2pnps(codonmat[, c(i, j)]))
                }
            }
            parallel::stopCluster(cl)
            OUT <- as.data.frame(OUT)
            OUT[["Ka"]] <- OUT[["dn"]]
            OUT[["Ks"]] <- OUT[["ds"]]
            OUT[["Ka/Ks"]] <- OUT[["dn/ds"]]
            attr(OUT, "model") <- "NG86"
            attr(OUT, "align") <- "FALSE"
            attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT)
        } else {
            if(.Platform$OS.type == "windows"){
                cl <- parallel::makeCluster(threads)
            }
            if(.Platform$OS.type != "windows"){
                cl <- parallel::makeForkCluster(threads)
            }
            doParallel::registerDoParallel(cl)
            i <- NULL
            j <- NULL
            aa <- local_cds2aa(cds, ...)
            cds <- local_cds2aa(cds, return.cds=TRUE, ...)
            OUT <- foreach(i = seq(from = 1, to = length(cds) - 1),
                .combine=rbind, .packages = c('foreach')) %dopar% {
                foreach(j = seq(from = i + 1, to = length(cds)),
                    .combine=rbind) %do% {
                    c(setNames(i, "Comp1"),
                    setNames(j, "Comp2"),
                    MSA2dist::codonmat2pnps(
                    MSA2dist::dnastring2codonmat(
                    local_cdsstring2codonaln(cds[c(i, j)], aa[c(i, j)],
                        ...))))
                }
            }
            parallel::stopCluster(cl)
            if(is.null(dim(OUT))){
                OUT <- as.data.frame(t(OUT))
            } else {
                OUT <- as.data.frame(OUT)
            }
            OUT[["Ka"]] <- OUT[["dn"]]
            OUT[["Ks"]] <- OUT[["ds"]]
            OUT[["Ka/Ks"]] <- OUT[["dn/ds"]]
            attr(OUT, "model") <- "NG86"
            attr(OUT, "align") <- "TRUE"
            attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT)
        }
    }
    if(model != "Li" && model != "NG86"){
        if(isMSA){
            OUT <- rcpp_KaKs(cdsstr = as.character(cds),
                sgc = sgc, method = model, verbose = verbose)
            OUT <- as.data.frame(t(
                tibble::column_to_rownames(tidyr::as_tibble(
                setNames(stringr::str_split(OUT$results_vec,
                "\t"), OUT$results_names)) %>%
                tibble::add_column(rownames=OUT$rownames),
                "rownames")))
            attr(OUT, "model") <- model
            attr(OUT, "align") <- "FALSE"
            attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT)
        } else{
            if(.Platform$OS.type == "windows"){
                cl <- parallel::makeCluster(threads)
            }
            if(.Platform$OS.type != "windows"){
                cl <- parallel::makeForkCluster(threads)
            }
            doParallel::registerDoParallel(cl)
            i <- NULL
            j <- NULL
            aa <- local_cds2aa(cds, ...)
            cds <- local_cds2aa(cds, return.cds=TRUE, ...)
            OUT <- foreach(i = seq(from = 1, to = length(cds) - 1),
                .combine=rbind, .packages = c('foreach')) %dopar% {
                foreach(j = seq(from = i + 1, to = length(cds)),
                    .combine=rbind) %do% {
                    tmp_out <- rcpp_KaKs(
                        cdsstr = as.character(
                        local_cdsstring2codonaln(cds[c(i, j)], aa[c(i, j)],
                            ...)),
                        sgc = sgc,
                        method = model)
                    tmp_out <- as.data.frame(t(
                        tibble::column_to_rownames(
                        tidyr::as_tibble(
                        setNames(stringr::str_split(
                        tmp_out$results_vec,
                        "\t"), tmp_out$results_names)) %>%
                        tibble::add_column(
                        rownames=tmp_out$rownames),
                        "rownames")))
                    tmp_out["Comp1"] <- i
                    tmp_out["Comp2"] <- j
                    tmp_out
                }
            }
            parallel::stopCluster(cl)
            OUT <- as.data.frame(OUT)
            attr(OUT, "model") <- model
            attr(OUT, "align") <- "TRUE"
            attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT)
        }
    }
}
