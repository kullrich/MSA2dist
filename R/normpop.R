#' @title normpop
#' @name normpop
#' @description This function normalizes population representation
#' @param pop [mandatory] population assignment:
#' \itemize{
#'     \item character vector (same length as seqnames)
#'     \item named character vector (mapping seqnames to pop)
#'     \item integer vector (same length as seqnames)
#'     \item named integer vector (mapping index to pop)
#'     \item list (see \code{addpop2dnastring()})
#' }
#' @param seqnames sequence names [mandatory]
#' @return list containing:
#' \itemize{
#'     \item pop_idx
#'     \item pop_name
#'     \item pop_levels
#'     \item pop_lookup
#'     \item poplist.names
#'     \item poplist.integer
#' }
#' @importFrom methods is
#' @importFrom methods slot
#' @seealso \code{\link[MSA2dist]{addpop2string}}
#' @examples
#' ## load example sequence data
#' data("iupac", package="MSA2dist")
#' ## create poplist
#' poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = grep("Mmd.IRA", names(iupac)),
#'     AFG = grep("Mmm.AFG", names(iupac)))
#' normpop(pop=poplist, seqnames=names(iupac))
#' ## mxixing index and names
#' poplist <- list(FRA = names(iupac)[grep("Mmd.FRA", names(iupac))],
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = names(iupac)[grep("Mmd.IRA", names(iupac))],
#'     AFG = grep("Mmm.AFG", names(iupac)))
#' normpop(pop=poplist, seqnames=names(iupac))
#' ## leaving out some sequences which will be assigned as "unassigned"
#' poplist <- list(FRA = names(iupac)[grep("Mmd.FRA", names(iupac))],
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = names(iupac)[grep("Mmd.IRA", names(iupac))])
#' normpop(pop=poplist, seqnames=names(iupac))
#' ## create character vector
#' popchar <- c(rep("FRA", 8),
#'     rep("GER", 8),
#'     rep("IRA", 8),
#'     rep("AFG", 6))
#' normpop(pop=popchar, seqnames=names(iupac))
#' ## create named character vector
#' popcharn <- c(rep("FRA", 8),
#'     rep("GER", 8),
#'     rep("IRA", 8),
#'     rep("AFG", 6))
#' names(popcharn) <- names(iupac)
#' normpop(pop=popcharn, seqnames=names(iupac))
#' ## create integer vector
#' popint <- c(rep(1L, 8),
#'     rep(2L, 8),
#'     rep(3L, 8),
#'     rep(4L, 6))
#' normpop(pop=popint, seqnames=names(iupac))
#' ## create named integer vector
#' popintn <- c(rep(1L, 8),
#'     rep(2L, 8),
#'     rep(3L, 8),
#'     rep(4L, 6))
#' names(popintn) <- names(iupac)
#' normpop(pop=popintn, seqnames=names(iupac))
#' @export normpop
#' @author Kristian K Ullrich

normpop <- function(pop, seqnames){
    stopifnot("seqnames must be provided" = !is.null(seqnames))
    stopifnot("pop cannot be NULL" = !is.null(pop))
    n <- length(seqnames)
    pop_name <- rep(NA_character_, n)
    if(is.character(pop) &&
        is.null(names(pop)) &&
        length(pop) == n){
        pop_name <- as.character(pop)
    } else if(is.character(pop) && !is.null(names(pop))){
        pop_name <- pop[seqnames]
        if(any(is.na(pop_name))){
            stop("Some sequence names missing in pop")
        }
        pop_name <- as.character(pop_name)
    } else if(is.integer(pop) &&
        is.null(names(pop)) &&
        length(pop) == n){
        lev <- sort(unique(pop))
        pop_name <- paste0("pop", match(pop, lev))
    } else if(is.integer(pop) &&
        !is.null(names(pop))){
        pop_name_tmp <- pop[seqnames]
        if(any(is.na(pop_name_tmp))){
            stop("Some sequence names missing in pop")
        }
        lev <- sort(unique(pop_name_tmp))
        pop_name <- paste0("pop", match(pop_name_tmp, lev))
    } else if(is.list(pop)){
        for(i in seq_along(pop)){
            members <- pop[[i]]
            if(is.character(members)){
                idx <- match(members, seqnames)
                if(any(is.na(idx))){
                    stop("Unknown sequence name(s) in population")
                }
            } else{
                idx <- as.integer(members)
                if(any(idx < 1 | idx > n)){
                    stop("Invalid index in population")
                }
            }
            pop_name[idx] <- names(pop)[i]
        }
    } else{
        stop("Unsupported pop format")
    }
    if(any(is.na(pop_name))){
        pop_name[is.na(pop_name)] <- "unassigned"
    }
    pop_levels <- unique(pop_name)
    pop_lookup <- setNames(pop_levels, seq_along(pop_levels))
    pop_idx <- match(pop_name, pop_levels)
    names(pop_idx) <- seqnames
    names(pop_name) <- seqnames
    poplist.names <- split(seqnames, pop_name)
    poplist.integer <- split(seq_len(n), pop_name)
    out <- list(
        pop_idx = pop_idx,
        pop_name = pop_name,
        pop_levels = pop_levels,
        pop_lookup = pop_lookup,
        poplist.names = poplist.names,
        poplist.integer = poplist.integer
    )
    return(out)
}
