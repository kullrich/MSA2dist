#' @title globalDeletion
#' @name globalDeletion
#' @description This function returns a \code{DNAStringSet} reduced by all
#' sites containing any gaps ("-", "+", ".") or missing ("N") sites.
#' @return \code{DNAStringSet}
#' @importFrom methods is slot
#' @importFrom Biostrings consensusMatrix
#' @param dna \code{DNAStringSet} [mandatory]
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' globalDeletion(cds1.cds2.aln)
#' @export globalDeletion
#' @author Kristian K Ullrich

globalDeletion<-function(dna){
    stopifnot("Error: input needs to be a DNAStringSet"=
        methods::is(dna, "DNAStringSet"))
    cM <- Biostrings::consensusMatrix(dna)
    globalDeletionSites <- which(apply(cM, 2, function(x) sum(x[15:18]) >= 1))
    if(length(globalDeletionSites) == 0){
        return(dna)
    }
    return(MSA2dist::dnabin2dnastring(
        MSA2dist::dnastring2dnabin(dna)[, -(globalDeletionSites)]))
}
