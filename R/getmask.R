#' @title getmask
#' @name getmask
#' @description This function shows the mask slot from a
#' \code{DNAStringSet} or an \code{AAStringSet}
#' \code{metadata} information.\cr
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @return \code{IRanges} information from \code{metadata}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @seealso \code{\link[MSA2dist]{addpop2string}}
#' @examples
#' ## load example sequence data
#' data(iupac, package="MSA2dist")
#' iupac.aa <- iupac |> cds2aa(shorten = TRUE)
#' ## create mask
#' mask1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
#' ## add mask
#' iupac.aa <- iupac.aa |> addmask2string(mask=mask1)
#' #(iupac.aa |> slot("metadata"))$mask
#' iupac.aa |> getmask()
#' @export getmask
#' @author Kristian K Ullrich

getmask <- function(seq){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    return(methods::slot(seq, "metadata")$mask)
}
