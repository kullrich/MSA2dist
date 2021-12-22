#' @title aastring2aabin
#' @name aastring2aabin
#' @description This function converts a \code{AAStringSet} into an \code{ape}
#' \code{DNAbin}.
#' @param aa \code{AAStringSet} [mandatory]
#' @return An object of class \code{DNAbin}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom ape as.AAbin
#' @seealso \code{\link[seqinr]{as.alignment}}
#' \code{\link[ape]{as.DNAbin.alignment}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' ## convert into AAbin
#' #aastring2aabin(cds2aa(cds1.cds2.aln))
#' cds1.cds2.aln |> cds2aa() |> aastring2aabin()
#' @export aastring2aabin
#' @author Kristian K Ullrich

aastring2aabin <- function(aa){
    stopifnot("Error: input needs to be an AAStringSet"=
                methods::is(aa, "AAStringSet"))
    return(ape::as.AAbin(aa))
}
