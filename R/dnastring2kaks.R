#' @title dnastring2kaks
#' @name dnastring2kaks
#' @description This function calculates Ka/Ks (pN/pS; according to
#' \emph{Li (1993)} or \emph{Nei and Gojobori (1986)} for all combinations of
#' a \code{DNAStringSet}.
#' @param cds \code{DNAStringSet} coding sequence alignment [mandatory]
#' @param model specify codon model either "Li" or "NG86" [default: Li]
#' @param threads number of parallel threads [default: 1]
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
#' @seealso \code{\link[seqinr]{kaks}}
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
#' #dnastring2kaks(hiv, model="Li")
#' hiv |> dnastring2kaks(model="Li")
#' #dnastring2kaks(hiv, model="NG86")
#' hiv |> dnastring2kaks(model="NG86")
#' #dnastring2kaks(hiv, model="NG86", threads=2)
#' hiv |> dnastring2kaks(model="NG86", threads=2)
#' @export dnastring2kaks
#' @author Kristian K Ullrich

dnastring2kaks <- function(cds,
    model = "Li",
    threads = 1){
    stopifnot("Error: input needs to be a DNAStringSet"=
        methods::is(cds, "DNAStringSet"))
    Comp1 <- FALSE
    Comp2 <- FALSE
    stopifnot("Error: either choose model 'Li' or 'NG86'"=
        model %in% c("Li", "NG86"))
    if(model == "Li"){
        if("Comp1" %in% names(cds)){
            names(cds)[which(names(cds) == "Comp1")] <- "_Comp1"
            Comp1 <- TRUE
        }
        if("Comp2" %in% names(cds)){
            names(cds)[which(names(cds) == "Comp2")] <- "_Comp2"
            Comp2 <- TRUE
        }
        OUT <- seqinr::kaks(dnastring2aln(cds))
        OUT.ka <- as.matrix(OUT$ka)
        colnames(OUT.ka) <- rownames(OUT.ka) <- gsub(" ", "", colnames(OUT.ka))
        OUT.ks <- as.matrix(OUT$ks)
        colnames(OUT.ks) <- rownames(OUT.ks) <- gsub(" ", "", colnames(OUT.ks))
        OUT.vka <- as.matrix(OUT$vka)
        colnames(OUT.vka) <- rownames(OUT.vka) <- gsub(" ", "",
            colnames(OUT.vka))
        OUT.vks <- as.matrix(OUT$vks)
        colnames(OUT.vks) <- rownames(OUT.vks) <- gsub(" ", "",
            colnames(OUT.vks))
        OUT.ka <- OUT.ka %>% tidyr::as_tibble(rownames = "Comp1") %>%
            tidyr::pivot_longer(colnames(OUT.ka),
            names_to = "Comp2", values_to = "ka") %>%
            dplyr::slice(uptriidx(ncol(OUT.ka)))
        OUT.ks <- OUT.ks %>% tidyr::as_tibble(rownames = "Comp1") %>%
            tidyr::pivot_longer(colnames(OUT.ks),
            names_to = "Comp2", values_to = "ks") %>%
            dplyr::slice(uptriidx(ncol(OUT.ks)))
        OUT.vka <- OUT.vka %>% tidyr::as_tibble(rownames = "Comp1") %>%
            tidyr::pivot_longer(colnames(OUT.vka),
            names_to = "Comp2", values_to = "vka") %>%
            dplyr::slice(uptriidx(ncol(OUT.vka)))
        OUT.vks <- OUT.vks %>% tidyr::as_tibble(rownames = "Comp1") %>%
            tidyr::pivot_longer(colnames(OUT.vks),
            names_to = "Comp2", values_to = "vks") %>%
            dplyr::slice(uptriidx(ncol(OUT.vks)))
        OUT <- OUT.ka %>% dplyr::left_join(OUT.ks) %>%
            dplyr::left_join(OUT.vka) %>% dplyr::left_join(OUT.vks)
        OUT <- as.data.frame(OUT)
        attr(OUT, "model") <- "Li"
        attr(OUT, "align") <- "FALSE"
        attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
        if(Comp1){
            OUT$Comp1 <- gsub("_Comp1", "Comp1", OUT$Comp1)
            OUT$Comp2 <- gsub("_Comp1", "Comp1", OUT$Comp2)
        }
        if(Comp2){
            OUT$Comp1 <- gsub("_Comp2", "Comp2", OUT$Comp1)
            OUT$Comp2 <- gsub("_Comp2", "Comp2", OUT$Comp2)
        }
        return(OUT)
    }
    if(model == "NG86"){
        #doMC::registerDoMC(threads)
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
                c(setNames(i, "Comp1"), setNames(j, "Comp2"),
                    MSA2dist::codonmat2pnps(codonmat[, c(i, j)]))
            }
        }
        parallel::stopCluster(cl)
        OUT <- as.data.frame(OUT)
        attr(OUT, "model") <- "NG86"
        attr(OUT, "align") <- "FALSE"
        attr(OUT, "MSA2dist.class") <- "dnastring2kaks"
        return(OUT)
    }
}
