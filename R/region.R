#' @title region
#' @name region
#' @description This function shows the region slot from a
#' \code{DNAStringSet} or an \code{AAStringSet}
#' \code{metadata} information.\cr
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @return region \code{IRanges} object from \code{metadata}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @seealso \code{\link[MSA2dist]{addpop2string}}
#' @examples
#' ## load example sequence data
#' data(iupac, package="MSA2dist")
#' iupac.aa <- iupac |> cds2aa(shorten = TRUE)
#' ## create region
#' region1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
#' ## add region
#' iupac.aa <- iupac.aa |> addregion2string(region=region1)
#' iupac.aa |> region()
#' @export region
#' @author Kristian K Ullrich

region <- function(seq){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    return(methods::slot(seq, "metadata")$region)
}
