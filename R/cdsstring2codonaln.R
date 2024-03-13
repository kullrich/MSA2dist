#' @title cdsstring2codonaln
#' @name cdsstring2codonaln
#' @description This function takes two sequences as \code{DNAStringSet},
#' and their corresponding \code{AAStringSet}, calculates
#' a global alignment and converts this alignment back into a codon alignment.
#' @param cds two sequences \code{DNAStringSet} [mandatory]
#' @param aa two sequences \code{AAStringSet} [mandatory]
#' @param type type of alignment (see
#' \code{\link[Biostrings]{pairwiseAlignment}}) [default: global]
#' @param substitutionMatrix substitution matrix representing the fixed
#' substitution scores for an alignment (see
#' \code{\link[Biostrings]{pairwiseAlignment}}) [default: BLOSUM62]
#' @param gapOpening the cost for opening a gap in the alignment (see
#' \code{\link[Biostrings]{pairwiseAlignment}}) [default: 10]
#' @param gapExtension the incremental cost incurred along the length of the
#' gap in the alignment (see \code{\link[Biostrings]{pairwiseAlignment}})
#' [default: 0.5]
#' @param remove.gaps specify if gaps in the codon alignment should be removed
#' [default: FALSE]
#' @return codon alignment as \code{DNAStringSet}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom methods is slot
#' @references Pagès, H et al. (2014) Biostrings: Efficient manipulation of
#' biological strings. \emph{R package version}, \bold{2(0)}.
#' @seealso \code{\link[Biostrings]{pairwiseAlignment}}
#' @examples
#' ## define two cds sequences
#' cds <- Biostrings::DNAStringSet(c("ATGCAACATTGC", "ATGCATTGC"))
#' names(cds) <- c("cds1", "cds2")
#' ## get protein alignment
#' aa <- MSA2dist::cds2aa(cds)
#' cdsstring2codonaln(cds, aa)
#' @export cdsstring2codonaln
#' @author Kristian K Ullrich

cdsstring2codonaln <- function(cds, aa, type="global",
    substitutionMatrix="BLOSUM62", gapOpening=10, gapExtension=0.5,
    remove.gaps=FALSE){
    stopifnot("Error: cds needs to be DNAStringSet"=
        {methods::is(cds, "DNAStringSet")})
    stopifnot("Error: aa needs to be either AAStringSet"=
        {methods::is(aa, "AAStringSet")})
    if(methods::is(cds, "DNAStringSet")){
        stopifnot("Error: cds needs to only contain two sequences"=
            length(cds) == 2)
    }
    if(methods::is(aa, "AAStringSet")){
        stopifnot("Error: aa needs to only contain two sequences"=
            length(aa) == 2)
    }
    xy.aln <- makePostalignedSeqs(Biostrings::pairwiseAlignment(aa[1], aa[2],
        type=type, substitutionMatrix=substitutionMatrix, gapOpening=gapOpening,
        gapExtension=gapExtension))[[1L]]
    names(xy.aln) <- names(aa)
    xy.cds.aln <- MSA2dist::pal2nal(xy.aln, cds, remove.gaps=remove.gaps)
    return(xy.cds.aln)
}
