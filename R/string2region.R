#' @title string2region
#' @name string2region
#' @description This function subsets a \code{DNAStringSet} or an
#' \code{AAStringSet} by a \code{mask} and \code{region} given one or both
#' options as \code{IRanges}.
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @param mask \code{IRanges} object indicating masked sites [default: NULL]
#' @param region \code{IRanges} object indicating region to use for dist
#' calculation (by default all sites are used) [default: NULL]
#' @param add indicate if mask and region should be added to \code{metadata}
#' [default: TRUE]
#' @return A \code{list} object with the following components:\cr
#' \code{DNAStringSet} or \code{AAStringSet}\cr
#' \code{regionUsed}\cr
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom IRanges IRanges IRangesList reduce start end findOverlaps
#' disjoin overlapsRanges
#' @seealso \code{\link[MSA2dist]{dnastring2dist}}
#' @examples
#' ## load example sequence data
#' data("hiv", package="MSA2dist")
#' ## create mask
#' mask1 <- IRanges::IRanges(start=c(11,41,71), end=c(20,50,80))
#' ## use mask
#' hiv.region <- hiv |> cds2aa() |> string2region(mask=mask1)
#' #(hiv.region |> slot("metadata"))$regionUsed
#' hiv.region |> regionused()
#' ## use region
#' region1 <- IRanges::IRanges(start=c(1,75), end=c(45,85))
#' hiv.region <- hiv |> cds2aa() |> string2region(region=region1)
#' #(hiv.region |> slot("metadata"))$regionUsed
#' hiv.region |> regionused()
#' ## use mask and region
#' hiv.region <- hiv |> cds2aa() |> string2region(mask=mask1, region=region1)
#' #(hiv.region |> slot("metadata"))$regionUsed
#' hiv.region |> regionused()
#' @export string2region
#' @author Kristian K Ullrich

string2region <- function(seq,
    mask = NULL,
    region = NULL,
    add = TRUE){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    region.seq <- IRanges::IRanges(start=1, end=unique(width(seq)))
    if(!is.null(mask)){
        mask <- mask[!(IRanges::start(mask)>unique(width(seq)))]
        IRanges::start(mask[IRanges::start(mask)<1]) <- 1
        IRanges::end(mask[IRanges::end(mask)>unique(width(seq))]) <-
            unique(width(seq))
        mask <- IRanges::reduce(mask)
        region.seq.mask.overlaps <- IRanges::findOverlaps(region.seq,
            IRanges::reduce(mask))
        if(length(region.seq.mask.overlaps) != 0){
            region.seq.mask.disjoin <- IRanges::disjoin(c(region.seq, mask))
            region.seq <- IRanges::reduce(region.seq.mask.disjoin[
                -(IRanges::findOverlaps(region.seq.mask.disjoin, mask)@from)])
        }
    }
    if(!is.null(region)){
        region <- region[!(IRanges::start(region)>unique(width(seq)))]
        IRanges::start(region[IRanges::start(region)<1]) <- 1
        IRanges::end(region[IRanges::end(region)>unique(width(seq))])<-
            unique(width(seq))
        region <- IRanges::reduce(region)
        region.seq <- IRanges::reduce(IRanges::overlapsRanges(
            region.seq, region))
        stopifnot("Error: specified region already masked or outside"=
            length(region.seq) != 0)
    }
    if(!is.null(mask) || !is.null(region)){
        seq.region <- MSA2dist::subString(
            seq, IRanges::start(region.seq), IRanges::end(region.seq))
        methods::slot(seq.region, "metadata") <- methods::slot(seq, "metadata")
        methods::slot(seq.region, "metadata")$regionUsed <- region.seq
        if(add){
            if(!is.null(mask)){
                MSA2dist::addmask2string(seq=seq.region, mask=mask,
                    append=FALSE)
            }
            if(!is.null(region)){
                MSA2dist::addregion2string(seq=seq.region, region=region,
                    append=FALSE)
            }
        }
    } else{seq.region <- seq
        methods::slot(seq.region, "metadata")$regionUsed <- region.seq
        if(!is.null(mask)){
            MSA2dist::addmask2string(seq=seq.region, mask=mask, append=FALSE)
        }
        if(!is.null(region)){
            MSA2dist::addregion2string(seq=seq.region, region=region,
                append=FALSE)
        }
    }
    return(seq.region)
}
