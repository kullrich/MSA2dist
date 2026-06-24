#' @title getpos
#' @name getpos
#' @description This function shows the position slot from a
#' \code{DNAStringSet} or an \code{AAStringSet}
#' \code{metadata} information.\cr
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @return \code{GenomicRanges} information from \code{metadata}
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
#' data(iupac, package="MSA2dist")
#' ## add position
#' iupac <- iupac |> addpos2string(chrom="chr1", start=1, end=1000)
#' #(iupac |> slot("metadata"))$GRanges
#' iupac |> getpos()
#' @export getpos
#' @author Kristian K Ullrich

getpos <- function(seq){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    return(methods::slot(seq, "metadata")$GRanges)
}
