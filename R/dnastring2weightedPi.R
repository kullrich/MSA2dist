#' @title dnastring2weightedPi
#' @name dnastring2weightedPi
#' @description This function calculates weightedPi for all
#' population combinations of a \code{DNAStringSet}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @param model specify model either "IUPAC" (diploid; decode IUPAC) or
#' "sequence" (haploid; only A,C,T,G counted as discrete alleles)
#' [default: IUPAC]
#' @param estimator specify estimator either
#' "count" (count-based; computes nucleotide diversity as the observed
#' proportion of pairwise allele differences, nDiff / nComp) or
#' "prob" (allele-frequency based; computes expected heterozygosity
#' as 1 - sum p_i^2) [default: count]
#' @param mask \code{IRanges} object indicating masked sites
#' [default: NULL]
#' @param region \code{IRanges} object indicating region to use for dist
#' calculation. Default is null, meaning all sites are used [default: NULL]
#' @param pop [default: NULL] optional population assignment:
#' \itemize{
#'     \item character vector (same length as dna)
#'     \item named character vector (mapping seqnames to pop)
#'     \item integer vector (same length as dna)
#'     \item named integer vector (mapping index to pop)
#'     \item list (see \code{addpop2dnastring()})
#' }
#' @return list containing:
#' \itemize{
#'     \item weightedPi
#'     \item populations
#'     \item model
#'     \item estimator
#'     \item regionUsed
#'     \item pop.info
#' }
#' @importFrom methods is
#' @importFrom methods slot
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom IRanges IRanges
#' @importFrom IRanges IRangesList
#' @importFrom IRanges reduce
#' @importFrom IRanges start
#' @importFrom IRanges end
#' @importFrom IRanges findOverlaps
#' @importFrom IRanges disjoin
#' @importFrom IRanges overlapsRanges
#' @seealso \code{\link[MSA2dist]{dnastring2dist}},
#' \code{\link[MSA2dist]{addpop2string}}
#' @examples
#' ## load example sequence data
#' data("iupac", package="MSA2dist")
#' #dnastring2weightedPi(iupac, model="IUPAC")
#' iupac |> dnastring2weightedPi(model="IUPAC")
#' #dnastring2weightedPi(iupac, model="IUPAC", estimator="prob")
#' iupac |> dnastring2weightedPi(model="IUPAC", estimator="prob")
#' #dnastring2weightedPi(iupac, model="sequence")
#' iupac |> dnastring2weightedPi(model="sequence")
#' #dnastring2weightedPi(iupac, model="sequence", estimator="prob")
#' iupac |> dnastring2weightedPi(model="sequence", estimator="prob")
#' ## create poplist
#' poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = grep("Mmd.IRA", names(iupac)),
#'     AFG = grep("Mmm.AFG", names(iupac)))
#' iupac <- iupac |> addpop2string(poplist)
#' #dnastring2weightedPi(iupac, pop=popinteger(iupac), model="IUPAC")
#' iupac |> dnastring2weightedPi(pop=popinteger(iupac), model="IUPAC")
#' #dnastring2weightedPi(iupac, pop=popnames(iupac), model="IUPAC")
#' iupac |> dnastring2weightedPi(pop=popnames(iupac), model="IUPAC")
#' #dnastring2weightedPi(iupac, pop=popinteger(iupac), model="IUPAC",
#' #estimator="prob")
#' iupac |> dnastring2weightedPi(pop=popinteger(iupac), model="IUPAC",
#' estimator="prob")
#' #dnastring2weightedPi(iupac, pop=popinteger(iupac), model="sequence")
#' iupac |> dnastring2weightedPi(pop=popinteger(iupac), model="sequence")
#' #dnastring2weightedPi(iupac, pop=popinteger(iupac), model="sequence",
#' #estimator="prob")
#' iupac |> dnastring2weightedPi(pop=popinteger(iupac), model="sequence",
#' estimator="prob")
#' ## create mask
#' mask1 <- IRanges::IRanges(start=c(1,61,121), end=c(30,90,150))
#' ## use mask
#' iupac |> dnastring2weightedPi(model="IUPAC", mask=mask1)
#' ## use region
#' region1 <- IRanges::IRanges(start=c(1,139), end=c(75,225))
#' iupac |> dnastring2weightedPi(model="IUPAC", region=region1)
#' ## use mask and region
#' iupac |> dnastring2weightedPi(model="IUPAC", mask=mask1, region=region1)
#' @export dnastring2weightedPi
#' @author Kristian K Ullrich

dnastring2weightedPi <- function(dna, model="IUPAC", estimator="count",
    mask=NULL, region=NULL, pop=NULL){
    stopifnot("Error: Input needs to be DNAStringSet"=
        methods::is(dna, "DNAStringSet"))
    if(!model %in% c("IUPAC", "sequence")){
        stop("model must be 'IUPAC' or 'sequence'")
    }
    if(!estimator %in% c("count", "prob")){
        stop("estimator must be 'count' or 'prob'")
    }
    region.dna <- IRanges::IRanges(start=1, end=unique(width(dna)))
    if(!is.null(mask) || !is.null(region)){
        dna.region <- MSA2dist::string2region(dna, mask=mask, region=region)
        dna.char <- as.character(dna.region)
        region.dna <- dna.region@metadata$regionUsed
    } else{dna.char <- as.character(dna)}
    dna.char.names <- make.unique(names(dna.char))
    n <- length(dna.char)
    if(is.null(pop)){
        pop.info <- MSA2dist::normpop(pop = rep("pop1", length(dna.char)),
            seqnames = dna.char.names)
        OUT <- rcpp_weightedPi(dnavector = dna.char, model = model,
            estimator = estimator)
        OUT$populations <- 1
    } else{
        pop.info <- MSA2dist::normpop(pop = pop, seqnames = dna.char.names)
        OUT <- rcpp_weightedPi_pop(dnavector = dna.char,
            pop_idx = pop.info$pop_idx, model = model, estimator = estimator)
    }
    OUT$regionUsed <- region.dna
    OUT$pop.info <- pop.info
    pop.names <- OUT$pop.info$pop_lookup[OUT$populations]
    names(OUT$weightedPi) <- pop.names
    names(OUT$populations) <- pop.names
    return(OUT)
}
