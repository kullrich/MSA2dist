#' @title addpop2string
#' @name addpop2string
#' @description This function adds population information to a
#' \code{DNAStringSet} or an \code{AAStringSet} and puts them into the
#' \code{metadata} information.\cr
#' __Note__: All unassigned sequences will be put into pop "unassigned"!\cr
#' Do not use "unassigned" as a population name!\cr
#' __Note__: Names in a population in the poplist must match sequence names!\cr
#' __Note__: Duplicated assignments are allowed!\cr
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @param poplist named \code{list} of populations either as index or names per
#' population (do not mix index and names in one population) [mandatory]
#' @return An object of class \code{DNAStringSet} or \code{AAStringSet}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @seealso \code{\link[MSA2dist]{addmask2string}},
#' \code{\link[MSA2dist]{addregion2string}},
#' \code{\link[MSA2dist]{addpos2string}}
#' @examples
#' ## load example sequence data
#' data(iupac, package="MSA2dist")
#' iupac.aa <- iupac |> cds2aa(shorten = TRUE)
#' ## create poplist
#' poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = grep("Mmd.IRA", names(iupac)),
#'     AFG = grep("Mmm.AFG", names(iupac)))
#' iupac.aa <- iupac.aa |> addpop2string(poplist)
#' #(iupac.aa |> slot("metadata"))$pop.integer
#' iupac.aa |> popinteger()
#' #(iupac.aa |> slot("metadata"))$pop.names
#' iupac.aa |> popnames()
#' ## mxixing index and names
#' poplist <- list(FRA = names(iupac)[grep("Mmd.FRA", names(iupac))],
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = names(iupac)[grep("Mmd.IRA", names(iupac))],
#'     AFG = grep("Mmm.AFG", names(iupac)))
#' iupac.aa <- iupac.aa |> addpop2string(poplist)
#' iupac.aa |> popinteger()
#' iupac.aa |> popnames()
#' ## leaving out some sequences which will be assigned as "unassigned"
#' poplist <- list(FRA = names(iupac)[grep("Mmd.FRA", names(iupac))],
#'     GER = grep("Mmd.GER", names(iupac)),
#'     IRA = names(iupac)[grep("Mmd.IRA", names(iupac))])
#' iupac.aa <- iupac.aa |> addpop2string(poplist)
#' iupac.aa |> popinteger()
#' iupac.aa |> popnames()
#' @export addpop2string
#' @author Kristian K Ullrich

addpop2string <- function(seq, poplist){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    stopifnot("Error: input needs to be a list"=
        methods::is(poplist, "list"))
    poplist.integer <- poplist
    poplist.names <- poplist
    poplist.type_is_integer <- which(unlist(lapply(poplist, is.integer)))
    poplist.type_is_names <- which(!(names(poplist)
        %in% names(poplist.type_is_integer)))
    poplist.integer[poplist.type_is_names] <- lapply(
        poplist[c(poplist.type_is_names)], function(x) unlist(lapply(x,
        function(y) grep(y, names(seq)))))
    poplist.names[poplist.type_is_integer] <- lapply(
        poplist[c(poplist.type_is_integer)], function(x) unlist(lapply(x,
        function(y) names(seq)[y])))
    unassigned.integer <- which(!(seq(from = 1, to = length(seq))
        %in% unlist(poplist.integer)))
    unassigned.names <- names(seq)[unassigned.integer]
    if(length(unassigned.names)!=0){
        poplist.integer$unassigned <- unassigned.integer
        poplist.names$unassigned <- unassigned.names
    }
    methods::slot(seq, "metadata")$pop.integer <- poplist.integer
    methods::slot(seq, "metadata")$pop.names <- poplist.names
    return(seq)
}
