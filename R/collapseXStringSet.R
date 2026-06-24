#' @title collapseXStringSet
#' @name collapseXStringSet
#' @description This function collapse identical sequences into
#' unique patterns with a mapping back to original order from
#' a \code{BStringSet}, \code{DNAStringSet}, \code{RNAStringSet}
#' or \code{AAStringSet}.
#' @param xstringset \code{XStringSet} [mandatory]
#' @return list containing:
#' \itemize{
#'     \item unique
#'     \item xmap
#' }
#' @importFrom methods is
#' @importFrom methods slot
#' @importFrom Biostrings BStringSet
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings RNAStringSet
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings readBStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readRNAStringSet
#' @importFrom Biostrings readAAStringSet
#' @seealso \code{\link[Biostrings]{XStringSet-class}}
#' @examples
#' data(woodmouse, package="ape")
#' #collapseXStringSet(dnabin2dnastring(woodmouse))
#' woodmouse |> dnabin2dnastring() |> collapseXStringSet()
#' @export collapseXStringSet
#' @author Kristian K Ullrich

collapseXStringSet <- function(xstringset){
    stopifnot("Error: input needs to be a XStringSet"=
        methods::is(xstringset, "XStringSet"))
    xchar <- as.character(xstringset)
    unique_idx <- !duplicated(xchar)
    unique_x <- xstringset[unique_idx]
    xmap <- match(xchar, unique_x)
    list(
        unique = unique_x,
        xmap = xmap
    )
}
