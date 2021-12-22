#' @title popnames
#' @name popnames
#' @description This function shows the population names slot from a
#' \code{DNAStringSet} or an \code{AAStringSet}
#' \code{metadata} information.\cr
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @return population names from \code{metadata}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @seealso \code{\link[MSA2dist]{addpop2string}}
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
#' popnames(iupac.aa)
#' @export popnames
#' @author Kristian K Ullrich

popnames <- function(seq){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    return(methods::slot(seq, "metadata")$pop.names)
}
