#' @title dnastring2dxy
#' @name dnastring2dxy
#' @description This function calculates dxy and fst for all
#' population combinations of a \code{DNAStringSet}.
#' @details
#' Four estimator/model combinations are available:
#'
#' \strong{1. model = "IUPAC", estimator = "prob"}
#'
#' IUPAC ambiguity codes are decoded into allele counts. For example,
#' R contributes one A and one G allele, while N contributes one A,
#' one C, one G and one T allele. Allele frequencies are estimated as
#'
#' \deqn{
#' p_i = \frac{n_i}{\sum_j n_j}
#' }{
#' p_i = n_i / sum(n_j)
#' }
#'
#' Within-population diversity is calculated as expected heterozygosity:
#'
#' \deqn{
#' \pi = 1 - \sum_i p_i^2
#' }{
#' pi = 1 - sum(p_i^2)
#' }
#'
#' Between-population divergence is calculated from allele frequencies:
#'
#' \deqn{
#' d_{xy} = 1 - \sum_i p_{x,i} p_{y,i}
#' }{
#' dxy = 1 - sum(px_i * py_i)
#' }
#'
#' Mean values are obtained by averaging across all valid sites.
#'
#' \strong{2. model = "IUPAC", estimator = "count"}
#'
#' IUPAC ambiguity codes are decoded into allele counts as above.
#' Nucleotide diversity is estimated directly from pairwise allele counts:
#'
#' \deqn{
#' \pi =
#' \frac{
#' AC + AG + AT + CG + CT + GT
#' }{
#' \binom{n}{2}
#' }
#' }{
#' pi = (AC+AG+AT+CG+CT+GT)/(n*(n-1)/2)
#' }
#'
#' where \eqn{A,C,G,T} denote allele counts at a site and
#' \eqn{n} is the total number of decoded alleles.
#'
#' Between-population divergence is estimated from observed pairwise
#' differences between populations:
#'
#' \deqn{
#' d_{xy} =
#' \frac{
#' A_x(C_y+G_y+T_y)+
#' C_x(A_y+G_y+T_y)+
#' G_x(A_y+C_y+T_y)+
#' T_x(A_y+C_y+G_y)
#' }{
#' n_x n_y
#' }
#' }{
#' dxy = differing allele pairs / (nx * ny)
#' }
#'
#' \strong{3. model = "sequence", estimator = "prob"}
#'
#' Only unambiguous nucleotides A, C, G and T are used. Ambiguous
#' symbols are ignored. Allele frequencies are estimated from the
#' observed nucleotides:
#'
#' \deqn{
#' p_i = \frac{n_i}{n}
#' }{
#' p_i = n_i / n
#' }
#'
#' Within-population diversity and between-population divergence are:
#'
#' \deqn{
#' \pi = 1 - \sum_i p_i^2
#' }{
#' pi = 1 - sum(p_i^2)
#' }
#'
#' \deqn{
#' d_{xy} = 1 - \sum_i p_{x,i} p_{y,i}
#' }{
#' dxy = 1 - sum(px_i * py_i)
#' }
#'
#' \strong{4. model = "sequence", estimator = "count"}
#'
#' Only unambiguous nucleotides A, C, G and T are counted as alleles.
#' Ambiguous symbols are ignored. Nucleotide diversity is estimated
#' from pairwise differences:
#'
#' \deqn{
#' \pi =
#' \frac{
#' AC + AG + AT + CG + CT + GT
#' }{
#' \binom{n}{2}
#' }
#' }{
#' pi = (AC+AG+AT+CG+CT+GT)/(n*(n-1)/2)
#' }
#'
#' Between-population divergence is estimated from pairwise
#' differences between populations:
#'
#' \deqn{
#' d_{xy} =
#' \frac{
#' A_x(C_y+G_y+T_y)+
#' C_x(A_y+G_y+T_y)+
#' G_x(A_y+C_y+T_y)+
#' T_x(A_y+C_y+G_y)
#' }{
#' n_x n_y
#' }
#' }{
#' dxy = differing allele pairs / (nx * ny)
#' }
#'
#' For all estimator/model combinations,
#'
#' \deqn{
#' F_{ST} =
#' \frac{d_{xy} -
#' \frac{\pi_x + \pi_y}{2}}
#' {d_{xy}}
#' }{
#' fst = (dxy - ((pi_x + pi_y)/2)) / dxy
#' }
#'
#' where \eqn{\pi_x} and \eqn{\pi_y} are the mean within-population
#' nucleotide diversities and \eqn{d_{xy}} is the mean between-population
#' nucleotide divergence.
#' @param dna \code{DNAStringSet} [mandatory]
#' @param model specify model either "IUPAC" (diploid; decode IUPAC) or
#' "sequence" (haploid; only A,C,T,G counted as discrete alleles)
#' [default: IUPAC]
#' @param estimator specify estimator either
#' "count" (count-based; computes nucleotide diversity as the observed
#' proportion of pairwise allele differences, nDiff / nComp) or
#' "prob" (allele-frequency based; computes expected heterozygosity
#' as 1 - sum p_i^2) [default: prob]
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
#'     \item dxy
#'     \item fst
#'     \item pi_within
#'     \item pi_between
#'     \item valid_sites
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
#' #dnastring2dxy(iupac, model="IUPAC")
#' iupac |> dnastring2dxy(model="IUPAC")
#' #dnastring2dxy(iupac, model="IUPAC", estimator="count")
#' iupac |> dnastring2dxy(model="IUPAC", estimator="count")
#' #dnastring2dxy(iupac, model="sequence")
#' iupac |> dnastring2dxy(model="sequence")
#' #dnastring2dxy(iupac, model="sequence", estimator="count")
#' iupac |> dnastring2dxy(model="sequence", estimator="count")
#' ## create poplist
#' poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = grep("Mmd.IRA", names(iupac)),
#'     AFG = grep("Mmm.AFG", names(iupac)))
#' iupac <- iupac |> addpop2string(poplist)
#' #dnastring2dxy(iupac, pop=popinteger(iupac), model="IUPAC")
#' iupac |> dnastring2dxy(pop=popinteger(iupac), model="IUPAC")
#' #dnastring2dxy(iupac, pop=popnames(iupac), model="IUPAC")
#' iupac |> dnastring2dxy(pop=popnames(iupac), model="IUPAC")
#' #dnastring2dxy(iupac, pop=popinteger(iupac), model="IUPAC",
#' #estimator="count")
#' iupac |> dnastring2dxy(pop=popinteger(iupac), model="IUPAC",
#' estimator="count")
#' #dnastring2dxy(iupac, pop=popinteger(iupac), model="sequence")
#' iupac |> dnastring2dxy(pop=popinteger(iupac), model="sequence")
#' #dnastring2dxy(iupac, pop=popinteger(iupac), model="sequence",
#' #estimator="count")
#' iupac |> dnastring2dxy(pop=popinteger(iupac), model="sequence",
#' estimator="count")
#' ## create mask
#' mask1 <- IRanges::IRanges(start=c(1,61,121), end=c(30,90,150))
#' ## use mask
#' iupac |> dnastring2dxy(model="IUPAC", mask=mask1)
#' ## use region
#' region1 <- IRanges::IRanges(start=c(1,139), end=c(75,225))
#' iupac |> dnastring2dxy(model="IUPAC", region=region1)
#' ## use mask and region
#' iupac |> dnastring2dxy(model="IUPAC", mask=mask1, region=region1)
#' @export dnastring2dxy
#' @author Kristian K Ullrich

dnastring2dxy <- function(dna, model="IUPAC", estimator="prob",
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
        OUT <- rcpp_dxy_fst_pop(dnavector = dna.char,
            pop_idx = pop.info$pop_idx, model = model, estimator = estimator)
        OUT$populations <- 1
    } else{
        pop.info <- MSA2dist::normpop(pop = pop, seqnames = dna.char.names)
        OUT <- rcpp_dxy_fst_pop(dnavector = dna.char,
            pop_idx = pop.info$pop_idx, model = model, estimator = estimator)
    }
    OUT$regionUsed <- region.dna
    OUT$pop.info <- pop.info
    pop.names <- OUT$pop.info$pop_lookup[OUT$populations]
    colnames(OUT$dxy) <- pop.names
    rownames(OUT$dxy) <- pop.names
    colnames(OUT$fst) <- pop.names
    rownames(OUT$fst) <- pop.names
    colnames(OUT$pi_within) <- pop.names
    rownames(OUT$pi_within) <- pop.names
    colnames(OUT$pi_between) <- pop.names
    rownames(OUT$pi_between) <- pop.names
    colnames(OUT$valid_sites) <- pop.names
    rownames(OUT$valid_sites) <- pop.names
    names(OUT$populations) <- pop.names
    return(OUT)
}
