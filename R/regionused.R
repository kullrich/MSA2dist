#' @title regionused
#' @name regionused
#' @description This function shows the region used slot from a
#' \code{DNAStringSet} or an \code{AAStringSet}
#' \code{metadata} information.\cr
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @return population names from \code{metadata}
#' @importFrom methods is
#' @importFrom methods slot
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAString
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readAAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings width
#' @importFrom Biostrings subseq
#' @seealso \code{\link[MSA2dist]{addpop2string}}
#' @examples
#' ## load example sequence data
#' data("hiv", package="MSA2dist")
#' ## create mask
#' mask1 <- IRanges::IRanges(start=c(11,41,71), end=c(20,50,80))
#' ## use mask
#' hiv.region <- hiv |> cds2aa() |> string2region(mask=mask1)
#' #(hiv.region |> slot("metadata"))$regionUsed
#' hiv.region |> regionused()
#' @export regionused
#' @author Kristian K Ullrich

regionused <- function(seq){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    return(methods::slot(seq, "metadata")$regionUsed)
}
