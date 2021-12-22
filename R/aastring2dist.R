#' @title aastring2dist
#' @name aastring2dist
#' @description This function calculates pairwise distances for all combinations
#' of a \code{AAStringSet}.
#' @param aa \code{AAStringSet} [mandatory]
#' @param threads number of parallel threads [default: 1]
#' @param score \code{score matrix} use a score matrix to calculate distances
#' [mandatory]
#' @param mask \code{IRanges} object indicating masked sites
#' [default: NULL]
#' @param region \code{IRanges} object indicating region to use for dist
#' calculation (by default all sites are used) [default: NULL]
#' @return A \code{data.frame} of pairwise distance values
#' \code{distSTRING}, sites used \code{sitesUsed} and region used
#' \code{regionUsed}
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
#' #aastring2dist(cds2aa(hiv), score=granthamMatrix())
#' hiv |> cds2aa() |> aastring2dist(score=granthamMatrix())
#' ## create mask
#' mask1 <- IRanges::IRanges(start=c(11,41,71), end=c(20,50,80))
#' ## use mask
#' hiv |> cds2aa() |> aastring2dist(score=granthamMatrix(), mask=mask1)
#' ## use region
#' region1 <- IRanges::IRanges(start=c(1,75), end=c(45,85))
#' hiv |> cds2aa() |> aastring2dist(score=granthamMatrix(), region=region1)
#' ## use mask and region
#' hiv |> cds2aa() |> aastring2dist(score=granthamMatrix(),
#'     mask=mask1, region=region1)
#' @export aastring2dist
#' @author Kristian K Ullrich

aastring2dist <- function(aa,
    threads = 1,
    score = NULL,
    mask = NULL,
    region = NULL){
    stopifnot("Error: input needs to be an AAStringSet"=
        methods::is(aa, "AAStringSet"))
    stopifnot("Error: set score matrix e.g 'granthamMatrix()'"= !is.null(score))
    region.aa <- IRanges::IRanges(start=1, end=unique(width(aa)))
    if(!is.null(mask) || !is.null(region)){
        aa.region <- MSA2dist::string2region(aa, mask=mask, region=region)
        aa.char <- as.character(aa.region)
        region.aa <- methods::slot(aa.region, "metadata")$regionUsed
    } else{aa.char <- as.character(aa)}
    OUT <- MSA2dist::rcpp_distSTRING(dnavector=aa.char,
        scoreMatrix=score, ncores=threads)
    OUT$distSTRING <- as.data.frame(OUT$distSTRING)
    OUT$sitesUsed <- as.data.frame(OUT$sitesUsed)
    OUT$regionUsed <- region.aa
    return(OUT)
}
