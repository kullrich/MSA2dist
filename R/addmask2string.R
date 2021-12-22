#' @title addmask2string
#' @name addmask2string
#' @description This function adds mask information as an \code{IRanges} object,
#' \code{START} and \code{END} information, to a
#' \code{DNAStringSet} or an \code{AAStringSet} and puts them into the
#' \code{metadata} information.
#' This information can be used to restrict the distance calculation to
#' specific regions of the \code{DNAStringSet} or the \code{AAStringSet}.
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @param mask \code{IRanges} object [mandatory]
#' @param append indicate if mask should be appended or overwritten
#' [default: TRUE]
#' @return An object of class \code{DNAStringSet} or \code{AAStringSet}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom IRanges IRanges IRangesList reduce start end findOverlaps
#' disjoin overlapsRanges
#' @seealso \code{\link[MSA2dist]{addregion2string}},
#' \code{\link[MSA2dist]{addpop2string}},
#' \code{\link[MSA2dist]{addpos2string}}
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
#' ## append mask
#' mask2 <- IRanges::IRanges(start=c(21), end=c(30))
#' iupac.aa <- iupac.aa |> addmask2string(mask=mask2)
#' #(iupac.aa |> slot("metadata"))$mask
#' iupac.aa |> getmask()
#' ## overwrite mask
#' iupac.aa <- iupac.aa |> addmask2string(mask=mask2, append=FALSE)
#' #(iupac.aa |> slot("metadata"))$mask
#' iupac.aa |> getmask()
#' ## reduce by mask
#' #iupac.aa.region <- iupac.aa |> string2region(mask=
#' #    (iupac.aa |> slot("metadata"))$mask)
#' iupac.aa.region <- iupac.aa |> string2region(mask=
#'     getmask(iupac.aa))
#' #iupac.aa.region |> slot("metadata")
#' iupac.aa.region |> getmask()
#' @export addmask2string
#' @author Kristian K Ullrich

addmask2string <- function(seq, mask = NULL, append = TRUE){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    stopifnot("Error: set mask"= !is.null(mask))
    stopifnot("Error: mask needs to be an IRanges object"=
        methods::is(mask, "IRanges"))
    if(append){
        if(is.null(slot(seq, "metadata")$mask)){
            #seq@metadata$mask <- IRanges::reduce(mask)
            methods::slot(seq, "metadata")$mask <- IRanges::reduce(mask)
        } else{
            methods::slot(seq, "metadata")$mask <- IRanges::reduce(
                c(methods::slot(seq, "metadata")$mask, mask))
        }
    } else{
        methods::slot(seq, "metadata")$mask <- IRanges::reduce(mask)
    }
    return(seq)
}
