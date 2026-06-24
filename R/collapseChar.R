#' @title collapseChar
#' @name collapseChar
#' @description This function collapse identical character
#' sequences into unique patterns with a mapping back to the
#' original order.
#' @param xchar \code{character vector} [mandatory]
#' @return list containing:
#' \itemize{
#'     \item unique
#'     \item xmap
#' }
#' @importFrom methods is
#' @importFrom methods slot
#' @examples
#' data(woodmouse, package="ape")
#' #collapseChar(as.character(dnabin2dnastring(woodmouse)))
#' woodmouse |> dnabin2dnastring() |> as.character() |>
#' collapseChar()
#' @export collapseChar
#' @author Kristian K Ullrich

collapseChar <- function(xchar){
    stopifnot("Error: input needs to be a character vector"=
        methods::is(xchar, "character"))
    unique_idx <- !duplicated(xchar)
    unique_x <- xchar[unique_idx]
    xmap <- match(xchar, unique_x)
    list(
        unique = unique_x,
        xmap = xmap
    )
}
