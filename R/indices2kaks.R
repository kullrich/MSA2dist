#' @title indices2kaks
#' @name indices2kaks
#' @description This function calculates Ka/Ks (pN/pS)
#' for all combinations given in an indices \code{list} of a
#' \code{DNAStringSet}.
#' If the sequences in the \code{DNAStringSet} are not a multiple-sequence
#' alignment, pairwise codon alignments can be calculated on the fly.
#' Models used and implemented according to
#' \emph{Li (1993)} (using \code{seqinr}) or
#' \emph{Nei and Gojobori (1986)} (own implementation) or models from
#' \code{KaKs_Calculator2} ported to \code{MSA2dist} with \code{Rcpp}.
#' @param cds \code{DNAStringSet} coding sequence alignment [mandatory]
#' @param indices \code{list} list of indices to calculate Ks/Ks [mandatory]
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
#' ## create indices
#' idx <- list(c(2, 3), c(5,7,9))
#' #indices2kaks(hiv, idx, model="Li")
#' hiv |> indices2kaks(idx, model="Li")
#' #indices2kaks(hiv, idx, model="NG86")
#' hiv |> indices2kaks(idx, model="NG86")
#' #indices2kaks(hiv, idx, model="NG86", threads=2)
#' hiv |> indices2kaks(idx, model="NG86", threads=2)
#'
#' ## define three unaligned cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATGCATTGC")
#' cds3 <- Biostrings::DNAString("ATGCAATGC")
#' cds_sequences <- Biostrings::DNAStringSet(list(cds1, cds2, cds3))
#' names(cds_sequences) <- c("cds1", "cds2", "cds3")
#' ## create indices
#' idx <- list(c(1, 2), c(1,3))
#' ## set isMSA to FALSE to automatically create pairwise codon alignments
#' #indices2kaks(cds_sequences, idx, model="Li", isMSA=FALSE)
#' cds_sequences |> indices2kaks(idx, model="Li", isMSA=FALSE)
#' @export indices2kaks
#' @author Kristian K Ullrich

indices2kaks <- function(cds,
    indices,
    model = "Li",
    threads = 1,
    isMSA = TRUE,
    sgc = "1",
    verbose = FALSE,
    ...){
    stopifnot("Error: input needs to be a DNAStringSet"=
        methods::is(cds, "DNAStringSet"))
    aa <- MSA2dist::cds2aa(cds, ...)
    Comp1 <- FALSE
    Comp2 <- FALSE
    seq1 <- FALSE
    seq2 <- FALSE
    cds.names <- names(cds)
    cds.names <- gsub(" ", "_",
        gsub(" $", "", gsub("\\s+", " ", cds.names)))
    names(cds) <- cds.names
    names(aa) <-  cds.names
    stopifnot("Error: either choose model 'Li' or 'NG86'
        or KaKs_Calculator2 model"=
        model %in% c("Li", "NG86",
            "NG", "LWL", "LPB", "MLWL", "MLPB", "GY", "YN",
            "MYN", "MS", "MA", "GNG", "GLWL", "GLPB", "GMLWL",
            "GMLPB", "GYN", "GMYN"))
    if(model == "Li"){
        if("Comp1" %in% names(cds)){
            names(cds)[which(names(cds) == "Comp1")] <- "_Comp1"
            names(aa)[which(names(aa) == "Comp1")] <- "_Comp1"
            Comp1 <- TRUE
        }
        if("Comp2" %in% names(cds)){
            names(cds)[which(names(cds) == "Comp2")] <- "_Comp2"
            names(aa)[which(names(aa) == "Comp2")] <- "_Comp2"
            Comp2 <- TRUE
        }
        if("seq1" %in% names(cds)){
            names(cds)[which(names(cds) == "seq1")] <- "_seq1"
            names(aa)[which(names(aa) == "seq1")] <- "_seq1"
            seq1 <- TRUE
        }
        if("seq2" %in% names(cds)){
            names(cds)[which(names(cds) == "seq2")] <- "_seq2"
            names(aa)[which(names(aa) == "seq2")] <- "_seq2"
            seq2 <- TRUE
        }
        if(isMSA){
            if(.Platform$OS.type == "windows"){
                cl <- parallel::makeCluster(threads)
            }
            if(.Platform$OS.type != "windows"){
                cl <- parallel::makeForkCluster(threads)
            }
            doParallel::registerDoParallel(cl)
            i <- NULL
            OUT_LIST <- foreach(i = seq(from = 1, to = length(indices)),
                .combine=rbind, .packages = c('foreach', 'tidyr')) %dopar% {
                OUT <- seqinr::kaks(dnastring2aln(cds[indices[[i]]]))
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
                OUT
            }
            parallel::stopCluster(cl)
            OUT_LIST <- as.data.frame(OUT_LIST)
            attr(OUT_LIST, "model") <- "Li"
            attr(OUT_LIST, "align") <- "TRUE"
            attr(OUT_LIST, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT_LIST)
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
            k <- NULL
            OUT_LIST <- foreach(i = seq(from = 1, to = length(indices)),
                .combine=rbind, .packages = c('foreach', 'tidyr')) %dopar% {
                foreach(j = seq(from = 1, to = length(indices[[i]]) - 1),
                .combine=rbind, .packages = c('foreach')) %do% {
                foreach(k = seq(from = j + 1, to = length(indices[[i]])),
                    .combine=rbind) %do% {
                    c(setNames(indices[[i]][j], "Comp1"),
                    setNames(indices[[i]][k], "Comp2"),
                    setNames(cds.names[indices[[i]][j]], "seq1"),
                    setNames(cds.names[indices[[i]][k]], "seq2"),
                    unlist(seqinr::kaks(MSA2dist::dnastring2aln(
                        MSA2dist::cdsstring2codonaln(cds[c(indices[[i]][j],
                            indices[[i]][k])], aa[c(indices[[i]][j],
                            indices[[i]][k])],
                            ...)))))
                    }
                }
            }
            parallel::stopCluster(cl)
            OUT_LIST <- as.data.frame(OUT_LIST)
            attr(OUT_LIST, "model") <- "Li"
            attr(OUT_LIST, "align") <- "TRUE"
            attr(OUT_LIST, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT_LIST)
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
            k <- NULL
            OUT_LIST <- foreach(i = seq(from = 1, to = length(indices)),
                .combine=rbind, .packages = c('foreach', 'tidyr')) %dopar% {
                codonmat_i <- MSA2dist::dnastring2codonmat(cds[indices[[i]]])
                foreach(j = seq(from = 1, to = ncol(codonmat_i) - 1),
                    .combine=rbind, .packages = c('foreach')) %do% {
                    foreach(k = seq(from = j + 1, to = ncol(codonmat_i)),
                        .combine=rbind) %do% {
                        c(setNames(indices[[i]][j], "Comp1"),
                        setNames(indices[[i]][k], "Comp2"),
                        MSA2dist::codonmat2pnps(codonmat_i[, c(j, k)]))
                    }
                }
            }
            parallel::stopCluster(cl)
            OUT_LIST <- as.data.frame(OUT_LIST)
            attr(OUT_LIST, "model") <- "NG86"
            attr(OUT_LIST, "align") <- "FALSE"
            attr(OUT_LIST, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT_LIST)
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
            k <- NULL
            OUT_LIST <- foreach(i = seq(from = 1, to = length(indices)),
                .combine=rbind, .packages = c('foreach', 'tidyr')) %dopar% {
                foreach(j = seq(from = 1, to = length(indices[[i]]) - 1),
                    .combine=rbind, .packages = c('foreach')) %do% {
                    foreach(k = seq(from = j + 1, to = length(indices[[i]])),
                        .combine=rbind) %do% {
                        c(setNames(indices[[i]][j], "Comp1"),
                        setNames(indices[[i]][k], "Comp2"),
                        MSA2dist::codonmat2pnps(
                        MSA2dist::dnastring2codonmat(
                        MSA2dist::cdsstring2codonaln(cds[c(indices[[i]][j],
                            indices[[i]][k])], aa[c(indices[[i]][j],
                            indices[[i]][k])],
                        ...))))
                    }
                }
            }
            parallel::stopCluster(cl)
            OUT_LIST <- as.data.frame(OUT_LIST)
            attr(OUT_LIST, "model") <- "NG86"
            attr(OUT_LIST, "align") <- "TRUE"
            attr(OUT_LIST, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT_LIST)
        }
    }
    if(model != "Li" && model != "NG86"){
        if(isMSA){
            if(.Platform$OS.type == "windows"){
                cl <- parallel::makeCluster(threads)
            }
            if(.Platform$OS.type != "windows"){
                cl <- parallel::makeForkCluster(threads)
            }
            doParallel::registerDoParallel(cl)
            i <- NULL
            OUT_LIST <- foreach(i = seq(from = 1, to = length(indices)),
                .combine=rbind, .packages = c('foreach', 'tidyr')) %dopar% {
                OUT <- rcpp_KaKs(cdsstr = as.character(cds[indices[[i]]]),
                    sgc = sgc, method = model, verbose = verbose)
                OUT <- as.data.frame(t(
                    tibble::column_to_rownames(tidyr::as_tibble(
                    setNames(stringr::str_split(OUT$results_vec,
                    "\t"), OUT$results_names)) %>%
                    tibble::add_column(rownames=OUT$rownames),
                    "rownames")))
                OUT["Comp1"] <- indices[[i]][as.numeric(unlist(OUT["Comp1"]))]
                OUT["Comp2"] <- indices[[i]][as.numeric(unlist(OUT["Comp2"]))]
                OUT
            }
            parallel::stopCluster(cl)
            attr(OUT_LIST, "model") <- model
            attr(OUT_LIST, "align") <- "FALSE"
            attr(OUT_LIST, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT_LIST)
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
            k <- NULL
            OUT_LIST <- foreach(i = seq(from = 1, to = length(indices)),
                .combine=rbind, .packages = c('foreach', 'tidyr')) %dopar% {
                foreach(j = seq(from = 1, to = length(indices[[i]]) - 1),
                    .combine=rbind, .packages = c('foreach')) %do% {
                    foreach(k = seq(from = j + 1, to = length(indices[[i]])),
                        .combine=rbind) %do% {
                    tmp_out <- rcpp_KaKs(
                        cdsstr = as.character(
                        MSA2dist::cdsstring2codonaln(cds[c(indices[[i]][j],
                            indices[[i]][k])], aa[c(indices[[i]][j],
                            indices[[i]][k])],
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
                    tmp_out["Comp1"] <- indices[[i]][j]
                    tmp_out["Comp2"] <- indices[[i]][k]
                    tmp_out
                    }
                }
            }
            parallel::stopCluster(cl)
            OUT_LIST <- as.data.frame(OUT_LIST)
            attr(OUT_LIST, "model") <- model
            attr(OUT_LIST, "align") <- "TRUE"
            attr(OUT_LIST, "MSA2dist.class") <- "dnastring2kaks"
            return(OUT_LIST)
        }
    }
}
