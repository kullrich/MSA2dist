#' @title makePostalignedSeqs
#' @name makePostalignedSeqs
#' @description This function is a fork from an internal function from
#' \code{Biostrings}
#' @param x x
#' @return get internal function makePostalignedSeqs
#' @seealso \code{\link[Biostrings]{pairwiseAlignment}},
#' \link[MSA2dist]{cds2codonaln}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATGCATTGC")
#' makePostalignedSeqs(Biostrings::pairwiseAlignment(
#'     cds2aa(Biostrings::DNAStringSet(cds1)),
#'     cds2aa(Biostrings::DNAStringSet(cds2))))
#' @export makePostalignedSeqs
#' @author Kristian K Ullrich

makePostalignedSeqs <- get('.makePostalignedSeqs',
    envir=asNamespace('Biostrings'), inherits=FALSE)
