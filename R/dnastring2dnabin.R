#' @title dnastring2dnabin
#' @name dnastring2dnabin
#' @description This function converts a \code{DNAStringSet} into an
#' \code{ape} \code{DNAbin}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @return An object of class \code{DNAbin}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom seqinr as.alignment
#' @importFrom ape as.DNAbin.alignment
#' @seealso \code{\link[seqinr]{as.alignment}}
#' \code{\link[ape]{as.DNAbin.alignment}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' ## convert into DNAbin
#' #dnastring2dnabin(cds1.cds2.aln)
#' cds1.cds2.aln |> dnastring2dnabin()
#' @export dnastring2dnabin
#' @author Kristian K Ullrich

dnastring2dnabin <- function(dna){
    stopifnot("Error: input needs to be a DNAStringSet"=
        methods::is(dna, "DNAStringSet"))
    alignment.nb <- length(dna)
    alignment.nam <- names(dna)
    alignment.seq <- tolower(as.character(dna))
    names(alignment.seq) <- NULL
    alignment.com <- NA
    alignment <- list(alignment.nb, alignment.nam, alignment.seq, alignment.com)
    names(alignment) <- c("nb", "nam", "seq", "com")
    attr(alignment, "class") <- "alignment"
    alignment.bin <- ape::as.DNAbin.alignment(alignment)
    return(alignment.bin)
}
