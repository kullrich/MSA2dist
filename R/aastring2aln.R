#' @title aastring2aln
#' @name aastring2aln
#' @description This function converts a \code{AAStringSet} into an
#' \code{seqinr} \code{alignment}.
#' @param aa \code{AAStringSet} [mandatory]
#' @return An object of class \code{alignment} which is a list with the
#' following components:\cr
#' \code{nb} the number of aligned sequences\cr
#' \code{nam} a vector of strings containing the names of the aligned
#' sequences\cr
#' \code{seq} a vector of strings containing the aligned sequences\cr
#' \code{com} a vector of strings containing the commentaries for each sequence
#' or \code{NA} if there are no comments
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom seqinr as.alignment
#' @seealso \code{\link[seqinr]{as.alignment}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' #aastring2aln(cds2aa(cds1.cds2.aln))
#' cds1.cds2.aln |> cds2aa() |> aastring2aln()
#' @export aastring2aln
#' @author Kristian K Ullrich

aastring2aln <- function(aa){
    stopifnot("Error: input needs to be an AAStringSet"=
        methods::is(aa, "AAStringSet"))
    alignment.nb <- length(aa)
    alignment.nam <- names(aa)
    alignment.seq <- tolower(as.character(aa))
    names(alignment.seq) <- NULL
    alignment.com <- NA
    alignment <- list(alignment.nb, alignment.nam, alignment.seq, alignment.com)
    names(alignment) <- c("nb", "nam", "seq", "com")
    attr(alignment, "class") <- "alignment"
    return(alignment)
}
