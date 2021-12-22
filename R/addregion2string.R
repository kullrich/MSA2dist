#' @title addregion2string
#' @name addregion2string
#' @description This function adds region information as an \code{IRanges}
#' object, \code{START} and \code{END} information, to a
#' \code{DNAStringSet} or an \code{AAStringSet} and puts them into the
#' \code{metadata} information.
#' This information can be used to restrict the distance calculation to
#' specific regions of the \code{DNAStringSet} or the \code{AAStringSet}.
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @param region \code{IRanges} object [mandatory]
#' @param append indicate if region should be appended or overwritten
#' [default: TRUE]
#' @return An object of class \code{DNAStringSet} or \code{AAStringSet}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom IRanges IRanges IRangesList reduce start end findOverlaps
#' disjoin overlapsRanges
#' @seealso \code{\link[MSA2dist]{addmask2string}},
#' \code{\link[MSA2dist]{addpop2string}},
#' \code{\link[MSA2dist]{addpos2string}}
#' @examples
#' ## load example sequence data
#' data(iupac, package="MSA2dist")
#' iupac.aa <- iupac |> cds2aa(shorten = TRUE)
#' ## create region
#' region1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
#' ## add region
#' iupac.aa <- iupac.aa |> addregion2string(region=region1)
#' #(iupac.aa |> slot("metadata"))$region
#' iupac.aa |> region()
#' ## append region
#' region2 <- IRanges::IRanges(start=c(21), end=c(30))
#' iupac.aa <- iupac.aa |> addregion2string(region=region2)
#' #(iupac.aa |> slot("metadata"))$region
#' iupac.aa |> region()
#' ## overwrite region
#' iupac.aa <- iupac.aa |> addregion2string(region=region2, append=FALSE)
#' #(iupac.aa |> slot("metadata"))$region
#' iupac.aa |> region()
#' ## reduce by region
#' #iupac.aa.region <- iupac.aa |> string2region(region=
#' #    (iupac.aa |> slot("metadata"))$region)
#' iupac.aa.region <- iupac.aa |> string2region(region=
#'     region(iupac.aa))
#' #iupac.aa.region |> slot("metadata")
#' iupac.aa.region |> region()
#' @export addregion2string
#' @author Kristian K Ullrich

addregion2string <- function(seq, region = NULL, append = TRUE){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    stopifnot("Error: set region"= !is.null(region))
    stopifnot("Error: region needs to be an IRanges object"=
        methods::is(region, "IRanges"))
    if(append){
        if(is.null(seq@metadata$region)){
            methods::slot(seq, "metadata")$region <- IRanges::reduce(region)
        } else{
            methods::slot(seq, "metadata")$region <- IRanges::reduce(
                c(methods::slot(seq, "metadata")$region, region))
        }
    } else{
        methods::slot(seq, "metadata")$region <- IRanges::reduce(region)
    }
    return(seq)
}
