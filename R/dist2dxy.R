#' @title dist2dxy
#' @name dist2dxy
#' @description This function takes the results of the
#' \code{dnastring2dist} function to calculate dxy and fst for all
#' population combinations and optional RND, Gmin, RNDmin and RNDmax
#' if one specifies \code{popx}, \code{popy} and \code{popout}.
#' @details
#' The function computes pairwise genetic differentiation statistics
#' from a precomputed sequence distance matrix.
#'
#' Let \eqn{d_{ij}} denote the per-site sequence distance between
#' sequences \eqn{i} and \eqn{j}, and let \eqn{S_{ij}} be the number
#' of aligned sites used.
#'
#' \bold{Within-population diversity}
#'
#' For population \eqn{x} with sequences \eqn{i in x}:
#'
#' \deqn{
#' D_x = \frac{1}{\binom{n_X}{2}} \sum_{i<j in X} d_{ij}
#' }{Dx = mean pairwise distance within population x}
#'
#' \bold{Between-population divergence (\eqn{D_{xy}})}
#'
#' For populations \eqn{x} and \eqn{y}:
#'
#' \deqn{
#' D_{xy} = \frac{1}{n_x n_y} \sum_{i in x} \sum_{j in y} d_{ij}
#' }{Dxy = mean pairwise distance between populations x and y}
#'
#' \bold{Fixation index (\eqn{F_{ST}})}
#'
#' Using the distance-based estimator:
#'
#' \deqn{
#' F_{ST} = \frac{D_{xy} - \frac{D_x + D_y}{2}}{D_{xy}}
#' }{FST = (Dxy - (Dx + Dy)/2) / Dxy}
#'
#' \bold{Minimum divergence (\eqn{D_{min}})}
#'
#' \deqn{
#' D_{min}(x,y) = \min_{i in x, j in y} d_{ij}
#' }{Dmin = minimum pairwise distance between populations}
#'
#' \bold{Outgroup-corrected divergence}
#'
#' Let \eqn{o} be an outgroup population:
#'
#' \deqn{
#' D(x,o) = \frac{1}{n_x n_o} \sum_{i in x} \sum_{k in o} d_{ik}
#' }
#' \deqn{
#' D(y,o) = \frac{1}{n_y n_o} \sum_{j in y} \sum_{k in o} d_{jk}
#' }
#' \deqn{
#' D_{out} = \frac{D(x,o) + D(y,o)}{2}
#' }
#'
#' \bold{RND statistics}
#'
#' \deqn{
#' RND = \frac{D_{xy}}{D_{out}}
#' }{RND = Dxy / Dout}
#'
#' \deqn{
#' G_{min} = \frac{D_{min}(x,y)}{D_{xy}}
#' }{Gmin = Dmin / Dxy}
#'
#' \deqn{
#' RND_{min} = \frac{D_{min}(x,y)}{D_{out}}
#' }{RNDmin = Dmin / Dout}
#'
#' \bold{Interpretation}
#' \itemize{
#'     \item Dxy: average sequence divergence between populations
#'     \item FST: relative differentiation scaled by
#' within-population diversity
#'     \item Dmin: closest haplotype distance between populations
#'     \item RND: divergence normalized by outgroup distance
#'     \item Gmin: proportion of minimum divergence relative to
#' mean divergence
#'     \item RNDmin: minimum divergence scaled by outgroup distance
#'     \item RNDmax: max divergence scaled by outgroup distance
#' }
#' @param d list returned by \code{dnastring2dist()}
#' @param pop [default: NULL] optional population assignment:
#' \itemize{
#'     \item character vector (same length as dna)
#'     \item named character vector (mapping seqnames to pop)
#'     \item integer vector (same length as dna)
#'     \item named integer vector (mapping index to pop)
#'     \item list (see \code{addpop2dnastring()})
#' }
#' @param popx [default: NULL] Character string specifying
#' population x in the rooted population topology ((x,y),o)
#' @param popy [default: NULL] Character string specifying
#' population y in the rooted population topology ((x,y),o)
#' @param popout [default: NULL] Character string specifying
#' outgroup population in the rooted population topology ((x,y),o)
#' @return list containing:
#' \itemize{
#'     \item dxy pairwise matrix
#'     \item fst pairwise matrix
#'     \item pi_within pop matrix
#'     \item pi_between 3d matrix
#'     \item valid_sites 
#'     \item populations
#'     \item popx
#'     \item popy
#'     \item popout
#'     \item Dx
#'     \item Dy
#'     \item Dxy
#'     \item FSTxy
#'     \item Dmin
#'     \item Dmax
#'     \item Dx_out
#'     \item Dy_out
#'     \item Dout
#'     \item Gmin
#'     \item RND
#'     \item RNDmin
#'     \item RNDmax
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
#' @importFrom stats as.dist
#' @importFrom stats median
#' @importFrom stats var
#' @seealso \code{\link[MSA2dist]{dnastring2dist}},
#' \code{\link[MSA2dist]{addpop2string}}
#' @examples
#' ## load example sequence data
#' data("iupac", package="MSA2dist")
#' #d <- dnastring2dist(iupac)
#' #dist2dxy(d)
#' iupac |> dnastring2dist() |> dist2dxy()
#' ## create poplist
#' poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = grep("Mmd.IRA", names(iupac)),
#'     AFG = grep("Mmm.AFG", names(iupac)))
#' iupac <- iupac |> addpop2string(poplist)
#' #d <- dnastring2dist(iupac)
#' #dist2dxy(d, pop=popinteger(iupac))
#' iupac |> dnastring2dist() |> dist2dxy(pop=popinteger(iupac))
#' #d <- dnastring2dist(iupac)
#' #dist2dxy(d, pop=popinteger(iupac),
#' #popx="FRA", popy="IRA", popout="AFG")
#' iupac |> dnastring2dist() |> dist2dxy(pop=popinteger(iupac),
#' popx="FRA", popy="IRA", popout="AFG")
#' @export dist2dxy
#' @author Kristian K Ullrich

dist2dxy <- function(d, pop=NULL,
    popx=NULL, popy=NULL, popout=NULL){
    stopifnot(is.list(d))
    stopifnot(!is.null(d$distSTRING))
    stopifnot(!is.null(d$sitesUsed))
    distmat <- as.matrix(d$distSTRING)
    sitesmat <- as.matrix(d$sitesUsed)
    dna.char.names <- rownames(distmat)
    if(is.null(dna.char.names)){
        stop("Distance matrix must contain row names")
    }
    if(is.null(pop)){
        pop.info <- MSA2dist::normpop(pop = rep("pop1",
            length(dna.char.names)), seqnames = dna.char.names)
    } else{
        pop.info <- MSA2dist::normpop(pop = pop, seqnames = dna.char.names)
    }
    pop_lookup <- pop.info$pop_lookup
    npops <- length(pop_lookup)
    pop_idx <- pop.info$pop_idx
    stats_names <- c("mean", "median", "min", "max", "var", "mean_sites")
    get_pair_stats <- function(ix, jx){
        if(identical(ix, jx)){
            vals <- as.vector(stats::as.dist(distmat[ix, jx]))
            sitesv <- as.vector(stats::as.dist(sitesmat[ix, jx]))
        } else{
            vals <- as.vector(distmat[ix, jx])
            sitesv <- as.vector(sitesmat[ix, jx])
        }
        return(
            c(mean = mean(vals, na.rm=TRUE),
            median = stats::median(vals, na.rm=TRUE),
            min = min(vals, na.rm=TRUE),
            max = max(vals, na.rm=TRUE),
            var = stats::var(vals, na.rm=TRUE),
            mean_sites = mean(sitesv, na.rm=TRUE))
        )
    }
    dxy <- matrix(NA, npops, npops)
    fst <- matrix(NA, npops, npops)
    valid_sites <- matrix(NA, npops, npops)
    rownames(dxy) <- unname(pop_lookup)
    rownames(fst) <- unname(pop_lookup)
    rownames(valid_sites) <- unname(pop_lookup)
    colnames(dxy) <- unname(pop_lookup)
    colnames(fst) <- unname(pop_lookup)
    colnames(valid_sites) <- unname(pop_lookup)
    pi_within <- matrix(NA, npops, 6)
    colnames(pi_within) <- stats_names
    rownames(pi_within) <- unname(pop_lookup)
    pi_between <- array(NA, dim = c(npops, npops, 6))
    dimnames(pi_between) <- list(unname(pop_lookup),
        unname(pop_lookup), stats_names)
    Dx <- NA
    Dy <- NA
    Dxy <- NA
    FSTxy <- NA
    Dmin <- NA
    Dmax <- NA
    Dx_out <- NA
    Dy_out <- NA
    Dout <- NA
    RND <- NA
    Gmin <- NA
    RNDmin <- NA
    RNDmax <- NA
    for(i in seq_along(pop_lookup)){
        pop_i <- pop_lookup[i]
        pop_i_idx <- which(pop_idx == names(pop_i))
        i_stats <- get_pair_stats(pop_i_idx, pop_i_idx)
        pi_within[i, ] <- i_stats
        if(i < npops){
            for(j in seq(from=i+1, to=npops)){
                pop_j <- pop_lookup[j]
                pop_j_idx <- which(pop_idx == names(pop_j))
                if(length(pop_i_idx) == 0 || length(pop_j_idx) == 0) next
                j_stats <- get_pair_stats(pop_j_idx, pop_j_idx)
                ij_stats <- get_pair_stats(pop_i_idx, pop_j_idx)
                pi_between[i, j, ] <- ij_stats
                dxy[i, j] <- ij_stats["mean"]
                dxy[j, i] <- ij_stats["mean"]
                dx <- i_stats["mean"]
                dy <- j_stats["mean"]
                fst[i, j] <- (dxy[i, j] - (dx + dy) / 2) / dxy[i, j]
                fst[j, i] <- (dxy[i, j] - (dx + dy) / 2) / dxy[i, j]
                valid_sites[i, j] <- ij_stats["mean_sites"]
                valid_sites[j, i] <- ij_stats["mean_sites"]
            }
        }
    }
    if(!is.null(popx) && !is.null(popy)){
        if(!popx %in% pop_lookup){
            stop("popx not found")
        }
        if(!popy %in% pop_lookup){
            stop("popy not found")
        }
        popx_ <- pop_lookup[pop_lookup == popx]
        popx_idx <- which(pop_idx == names(popx_))
        popy_ <- pop_lookup[pop_lookup == popy]
        popy_idx <- which(pop_idx == names(popy_))
        Dx <- pi_within[popx, ]["mean"]
        Dy <- pi_within[popy, ]["mean"]
        Dxy <- pi_between[popx, popy, "mean"]
        Dmin <- pi_between[popx, popy, "min"]
        Dmax <- pi_between[popx, popy, "max"]
        FSTxy <- (Dxy - (Dx + Dy) / 2) / Dxy
        Gmin <- Dmin / Dxy
        if(!is.null(popout)){
            if(!popout %in% pop_lookup){
                stop("popout not found")
            }
            Dx_out <- dxy[popx, popout]
            Dy_out <- dxy[popy, popout]
            Dout <- (Dx_out + Dy_out) / 2
            RND <- Dxy / Dout
            RNDmin <- Dmin / Dout
            RNDmax <- Dmax /Dout
        }
    }
    OUT <- list(
        dxy = dxy,
        fst = fst,
        pi_within = pi_within,
        pi_between = pi_between,
        valid_sites = valid_sites,
        populations = setNames(as.numeric(names(pop_lookup)), pop_lookup),
        popx = popx,
        popy = popy,
        popout = popout,
        Dx = Dx,
        Dy = Dy,
        Dxy = Dxy,
        FSTxy = FSTxy,
        Dmin = Dmin,
        Dmax = Dmax,
        Dx_out =Dx_out,
        Dy_out = Dy_out,
        Dout = Dout,
        Gmin = Gmin,
        RND = RND,
        RNDmin = RNDmin,
        RNDmax = RNDmax
    )
    OUT$regionUsed <- d$regionUsed
    OUT$pop.info <- pop.info
    return(OUT)
}
