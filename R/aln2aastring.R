#' @title aln2aastring
#' @name aln2aastring
#' @description This function converts a \code{seqinr} \code{alignment} into
#' an \code{AAStringSet}.
#' @param aln \code{seqinr} \code{alignment} [mandatory]
#' @return An object of class \code{AAStringSet}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom seqinr as.alignment
#' @seealso \code{\link[seqinr]{as.alignment}}
#' \code{\link[Biostrings]{AAStringSet}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' #aastring2aln(cds2aa(cds1.cds2.aln))
#' cds1.cds2.aln |> cds2aa() |> aastring2aln() |> aln2aastring()
#' @export aln2aastring
#' @author Kristian K Ullrich

aln2aastring <- function(aln){
    aa <- Biostrings::AAStringSet(toupper(unlist(aln$seq)))
    names(aa) <- aln$nam
    return(aa)
}
