#' @title globalDeletionAA
#' @name globalDeletionAA
#' @description This function returns an \code{AAStringSet} reduced by all
#' sites containing any gaps ("-", "+", ".") or missing ("X") sites.
#' @return \code{AAStringSet}
#' @importFrom methods is
#' @importFrom methods slot
#' @importFrom Biostrings consensusMatrix
#' @param aa \code{AAStringSet} [mandatory]
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' #globalDeletionAA(cds2aa(cds1.cds2.aln))
#' cds1.cds2.aln |> cds2aa() |> globalDeletionAA()
#' @export globalDeletionAA
#' @author Kristian K Ullrich

globalDeletionAA<-function(aa){
    stopifnot("Error: input needs to be an AAStringSet"=
        methods::is(aa, "AAStringSet"))
    cM <- Biostrings::consensusMatrix(aa)
    globalDeletionSites <- which(apply(cM, 2, function(x) {
        sum(x[c(26, 28:30)])} >= 1))
    if(length(globalDeletionSites) == 0){
        return(aa)
    }
    return(MSA2dist::aabin2aastring(
        MSA2dist::aastring2aabin(aa)[, -(globalDeletionSites)]))
}
