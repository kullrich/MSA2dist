#' @title addpos2string
#' @name addpos2string
#' @description This function adds \code{GenomicRanges} information,
#' \code{CHROM}, \code{START} and \code{END} to a
#' \code{DNAStringSet} or an \code{AAStringSet} and puts them into the
#' \code{metadata} information.
#' This information can be used to find overlaps with a chromosome wide mask.
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @param chrom chromosome name [mandatory]
#' @param start start position [mandatory]
#' @param end end position [mandatory]
#' @return An object of class \code{DNAStringSet} or \code{AAStringSet}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges IRangesList reduce start end findOverlaps
#' disjoin overlapsRanges
#' @seealso \code{\link[MSA2dist]{addmask2string}},
#' \code{\link[MSA2dist]{addregion2string}},
#' \code{\link[MSA2dist]{addpop2string}}
#' @examples
#' ## load example sequence data
#' data(iupac, package="MSA2dist")
#' ## add position
#' iupac <- iupac |> addpos2string(chrom="chr1", start=1, end=1000)
#' #(iupac |> slot("metadata"))$GRanges
#' iupac |> getpos()
#' @export addpos2string
#' @author Kristian K Ullrich

addpos2string <- function(seq, chrom = NULL, start = NULL, end = NULL){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    stopifnot("Error: set chromosome name"= !is.null(chrom))
    stopifnot("Error: set start"= !is.null(start))
    stopifnot("Error: set end"= !is.null(end))
    if((end - start + 1)!=unique(width(seq))){
        warning("defined region (start, end) has unequal length to seq")
    }
    seq.GRanges <- GenomicRanges::GRanges(seqnames=chrom,
        ranges=IRanges::IRanges(start=start, end=end))
    methods::slot(seq, "metadata")$GRanges <- seq.GRanges
    return(seq)
}
