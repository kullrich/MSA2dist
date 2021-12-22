#' @title compareCodons
#' @name compareCodons
#' @description This function compares two codons and returns the number of syn
#'and non-syn sites according to \emph{Nei and Gojobori (1986)}.
#' @param codA \code{codon} A [mandatory]
#' @param codB \code{codon} B [mandatory]
#' @return vector of syn and non-syn sites
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom seqinr kaks
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
#' compareCodons("AAA","TTA")
#' compareCodons("AAA","TAT")
#' compareCodons("AAA","ATT")
#' compareCodons("AAA","TTT")
#' ## load example sequence data
#' data("hiv", package="MSA2dist")
#' compareCodons(dnastring2codonmat(hiv)[1,1], dnastring2codonmat(hiv)[1,2])
#' @export compareCodons
#' @author Kristian K Ullrich

compareCodons <- function(codA, codB){
    type <- "syn"
    if(codA == codB){
        return(c(0, 0, 0))
    }
    codAsplit <- strsplit(codA, "")[[1]]
    codBsplit <- strsplit(codB, "")[[1]]
    if("N" %in% codAsplit | "N" %in% codBsplit){
        type <- "Ns"
        return(c(NA, NA, NA))
    }
    if("-" %in% codAsplit | "-" %in% codBsplit){
        type <- "indel"
        return(c(0, 0, 1))
    }
    if(MSA2dist::GENETIC_CODE_TCAG[codA, 2] !=
        MSA2dist::GENETIC_CODE_TCAG[codB, 2]){
        type <- "nonsyn"
    }
    codAcodB_diff_idx <- which(codAsplit != codBsplit)
    if(length(codAcodB_diff_idx) == 1){
        if(type == "syn"){
            return(c(1, 0, 0))
        }
        if(type == "nonsyn"){
            return(c(0, 1, 0))
        }
    }
    if(length(codAcodB_diff_idx) == 2){
        ## AAA (Lys) -> TTA (Leu)
        ## AGC (Ser) -> TCC (Ser)
        if(all(codAcodB_diff_idx == c(1, 2))){
            codC <- paste0(codBsplit[1], codAsplit[2], codAsplit[3])
            codD <- paste0(codAsplit[1], codBsplit[2], codAsplit[3])
        }
        ## AAA -> TAT
        if(all(codAcodB_diff_idx == c(1, 3))){
            codC <- paste0(codBsplit[1], codAsplit[2], codAsplit[3])
            codD <- paste0(codAsplit[1], codAsplit[2], codBsplit[3])
        }
        ## AAA -> ATT
        if(all(codAcodB_diff_idx == c(2, 3))){
            codC <- paste0(codAsplit[1], codBsplit[2], codAsplit[3])
            codD <- paste0(codAsplit[1], codAsplit[2], codBsplit[3])
        }
        tmp_syn <- 0
        if(MSA2dist::GENETIC_CODE_TCAG[codA, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codC, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codA, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codD, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codB, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codC, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codB, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codD, 2]){
            tmp_syn <- tmp_syn + 1
        }
        tmp_syn <- tmp_syn / 2
        return(c(tmp_syn, 2 - tmp_syn, 0))
    }
    if(length(codAcodB_diff_idx) == 3){
        ## AAA -> TTT
        codC <- paste0(codBsplit[1], codAsplit[2], codAsplit[3])
        codD <- paste0(codAsplit[1], codBsplit[2], codAsplit[3])
        codE <- paste0(codAsplit[1], codAsplit[2], codBsplit[3])
        codF <- paste0(codBsplit[1], codBsplit[2], codAsplit[3])
        codG <- paste0(codBsplit[1], codAsplit[2], codBsplit[3])
        codH <- paste0(codAsplit[1], codBsplit[2], codBsplit[3])
        tmp_syn <- 0
        if(MSA2dist::GENETIC_CODE_TCAG[codA, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codC, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codA, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codD, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codA, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codE, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codB, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codF, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codB, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codG, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codB, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codH, 2]){
            tmp_syn <- tmp_syn + 1
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codC, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codF, 2]){
            tmp_syn <- tmp_syn + 0.5
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codC, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codG, 2]){
            tmp_syn <- tmp_syn + 0.5
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codD, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codF, 2]){
            tmp_syn <- tmp_syn + 0.5
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codD, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codH, 2]){
            tmp_syn <- tmp_syn + 0.5
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codE, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codG, 2]){
            tmp_syn <- tmp_syn + 0.5
        }
        if(MSA2dist::GENETIC_CODE_TCAG[codE, 2] ==
            MSA2dist::GENETIC_CODE_TCAG[codH, 2]){
            tmp_syn <- tmp_syn + 0.5
        }
        tmp_syn <- tmp_syn / 3
        return(c(tmp_syn, 3 - tmp_syn, 0))
    }
}
