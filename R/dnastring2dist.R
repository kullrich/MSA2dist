#' @title dnastring2dist
#' @name dnastring2dist
#' @description This function calculates pairwise distances for all
#' combinations of a \code{DNAStringSet}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @param model specify model either "IUPAC" or any model from
#' \code{ape::dist.dna} [default: IUPAC]
#' @param threads number of parallel threads [default: 1]
#' @param score \code{score matrix} use score matrix to calculate
#' distances [default: NULL]
#' @param mask \code{IRanges} object indicating masked sites
#' [default: NULL]
#' @param region \code{IRanges} object indicating region to use for dist
#' calculation. Default is null, meaning all sites are used [default: NULL]
#' @param ... other \code{ape::dist.dna} parameters
#' (see \code{\link[ape]{dist.dna}})
#' @return A data.frame of pairwise distance values \code{distSTRING} and
#' sites used \code{sitesUsed}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom ape dist.dna
#' @importFrom IRanges IRanges IRangesList reduce start end findOverlaps
#' disjoin overlapsRanges
#' @seealso \code{\link[ape]{dist.dna}}
#' @examples
#' ## load example sequence data
#' data("hiv", package="MSA2dist")
#' #dnastring2dist(hiv, model="IUPAC")
#' hiv |> dnastring2dist(model="IUPAC")
#' #dnastring2dist(hiv, model="K80")
#' hiv |> dnastring2dist(model="K80")
#' data("woodmouse", package="ape")
#' #dnastring2dist(dnabin2dnastring(woodmouse), score=iupacMatrix())
#' woodmouse |> dnabin2dnastring() |> dnastring2dist()
#' #dnastring2dist(hiv, model = "IUPAC", threads = 2)
#' hiv |> dnastring2dist(model = "IUPAC", threads = 2)
#' ## create mask
#' mask1 <- IRanges::IRanges(start=c(1,61,121), end=c(30,90,150))
#' ## use mask
#' hiv |> dnastring2dist(model="IUPAC", mask=mask1)
#' ## use region
#' region1 <- IRanges::IRanges(start=c(1,139), end=c(75,225))
#' hiv |> dnastring2dist(model="IUPAC", region=region1)
#' ## use mask and region
#' hiv |> dnastring2dist(model="IUPAC", mask=mask1, region=region1)
#' @export dnastring2dist
#' @author Kristian K Ullrich

dnastring2dist <- function(dna,
    model = "IUPAC",
    threads = 1,
    score = NULL,
    mask = NULL,
    region = NULL,
    ...){
    stopifnot("Error: Input needs to be DNAStringSet"=
        methods::is(dna, "DNAStringSet"))
    region.dna <- IRanges::IRanges(start=1, end=unique(width(dna)))
    if(!is.null(mask) || !is.null(region)){
        dna.region <- MSA2dist::string2region(dna, mask=mask, region=region)
        dna.char <- as.character(dna.region)
        region.dna <- dna.region@metadata$regionUsed
    } else{dna.char <- as.character(dna)}
    if(!is.null(score)){
        OUT <- rcpp_distSTRING(dnavector = dna.char,
            scoreMatrix = score, ncores = threads)
        OUT$distSTRING <- as.data.frame(OUT$distSTRING)
        OUT$sitesUsed <- as.data.frame(OUT$sitesUsed)
    }
    if(is.null(score)){
        stopifnot("Error: either choose model 'IUPAC' or '?ape::dist.dna'"=
            model %in%
            c("IUPAC", "raw", "N", "TS", "TV", "JC69", "K80", "F81",
            "K81","F84", "BH87", "T92", "TN93", "GG95", "logdet",
            "paralin", "indel", "indelblock"))
        if(model == "IUPAC"){
            OUT <- rcpp_distSTRING(dnavector = dna.char,
                scoreMatrix=iupacMatrix(), ncores=threads)
            OUT$distSTRING <- as.data.frame(OUT$distSTRING)
            OUT$sitesUsed <- as.data.frame(OUT$sitesUsed)
        }
        if(model != "IUPAC"){
            if(!is.null(mask) || !is.null(region)){
                distSTRING_ <- as.matrix(ape::dist.dna(
                    x=MSA2dist::dnastring2dnabin(dna.region$seq),
                    model=model, ...))
                sitesUsed_ <- rcpp_pairwiseDeletionDNA(
                    dnavector=dna.char,
                    ncores=threads)
            } else{
                distSTRING_ <- as.matrix(ape::dist.dna(
                    x=MSA2dist::dnastring2dnabin(dna),
                    model=model, ...))
                sitesUsed_ <- rcpp_pairwiseDeletionDNA(
                    dnavector=as.character(dna),
                    ncores=threads)
            }
            OUT <- list(distSTRING = as.data.frame(distSTRING_))
            OUT <- append(OUT, sitesUsed_)
            OUT$sitesUsed <- as.data.frame(OUT$sitesUsed)
        }
    }
    OUT$regionUsed <- region.dna
    return(OUT)
}
