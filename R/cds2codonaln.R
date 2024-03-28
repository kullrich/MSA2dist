#' @title cds2codonaln
#' @name cds2codonaln
#' @description This function takes two single sequence \code{DNAString}'s or
#' two single sequence \code{DNAStringSet}'s, converts them into aa, calculates
#' a global alignment and converts this alignment back into a codon alignment.
#' @param cds1 single sequence \code{DNAStringSet} or \code{DNAString}
#' [mandatory]
#' @param cds2 single sequence \code{DNAStringSet} or \code{DNAString}
#' [mandatory]
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
#' @param ... other cds2aa parameters
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
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATGCATTGC")
#' cds2codonaln(cds1, cds2)
#' @export cds2codonaln
#' @author Kristian K Ullrich

cds2codonaln <- function(cds1, cds2, type="global",
    substitutionMatrix="BLOSUM62", gapOpening=10, gapExtension=0.5,
    remove.gaps=FALSE, ...){
    stopifnot("Error: cds1 needs to be either DNAString or DNAStringSet"=
        {methods::is(cds1, "DNAString") || methods::is(cds1, "DNAStringSet")})
    stopifnot("Error: cds2 needs to be either DNAString or DNAStringSet"=
        {methods::is(cds2, "DNAString") || methods::is(cds2, "DNAStringSet")})
    if(methods::is(cds1, "DNAString")){
        x.aa <- MSA2dist::cds2aa(Biostrings::DNAStringSet(cds1), ...)[[1]]
        x.name <- "cds1"
        cds1 <- MSA2dist::cds2aa(Biostrings::DNAStringSet(cds1),
            return.cds=TRUE, ...)[[1]]
    }
    if(methods::is(cds2, "DNAString")){
        y.aa <- MSA2dist::cds2aa(Biostrings::DNAStringSet(cds2), ...)[[1]]
        y.name <- "cds2"
        cds2 <- MSA2dist::cds2aa(Biostrings::DNAStringSet(cds2),
            return.cds=TRUE, ...)[[1]]
    }
    if(methods::is(cds1, "DNAStringSet")){
        stopifnot("Error: cds1 needs to only contain one sequence"=
            length(cds1) == 1)
        x.aa <- MSA2dist::cds2aa(cds1, ...)[[1]]
        x.name <- names(cds1)
        cds1 <- MSA2dist::cds2aa(cds1, return.cds=TRUE, ...)[[1]]
    }
    if(methods::is(cds2, "DNAStringSet")){
        stopifnot("Error: cds2 needs to only contain one sequence"=
            length(cds2) == 1)
        y.aa <- MSA2dist::cds2aa(cds2, ...)[[1]]
        y.name <- names(cds2)
        cds2 <- MSA2dist::cds2aa(cds2, return.cds=TRUE, ...)[[1]]
    }
    xy.aln <- makePostalignedSeqs(Biostrings::pairwiseAlignment(x.aa, y.aa,
        type=type, substitutionMatrix=substitutionMatrix, gapOpening=gapOpening,
        gapExtension=gapExtension))[[1L]]
    names(xy.aln) <- c(x.name, y.name)
    if(type=="local"){
        xy.aln.x <- gsub("\\*", "X", gsub("-", "", xy.aln[[1]]))
        cds1_local_pos <- gregexpr(xy.aln.x, gsub("\\*", "X", x.aa))
        xy.aln.y <- gsub("\\*", "X", gsub("-", "", xy.aln[[2]]))
        cds2_local_pos <- gregexpr(xy.aln.y, gsub("\\*", "X", y.aa))
        cds1_local <- Biostrings::subseq(cds1,
            (cds1_local_pos[[1]][1]*3)-2,
            (cds1_local_pos[[1]][1]+nchar(xy.aln.x)-1)*3)
        cds2_local <- Biostrings::subseq(cds2,
            (cds2_local_pos[[1]][1]*3)-2,
            (cds2_local_pos[[1]][1]+nchar(xy.aln.y)-1)*3)
        xy.cds <- setNames(Biostrings::DNAStringSet(list(cds1_local,
            cds2_local)), c(x.name, y.name))
    } else{
        xy.cds <- setNames(Biostrings::DNAStringSet(list(cds1, cds2)),
            c(x.name, y.name))
    }
    xy.cds.aln <- MSA2dist::pal2nal(xy.aln, xy.cds, remove.gaps=remove.gaps)
    return(xy.cds.aln)
}
