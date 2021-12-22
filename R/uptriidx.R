#' @title uptriidx
#' @name uptriidx
#' @description This function returns upper tri index for usage with
#' \code{pivot_long} reduction.
#' @param n dimension of initial matrix [mandatory]
#' @param diag indicate if diag should be retained [default: FALSE]
#' @return list of positions
#' @examples
#' uptriidx(10)
#' @export uptriidx
#' @author Kristian K Ullrich

uptriidx <- function(n, diag = FALSE){
    stopifnot("Error: n needs to larger than 1"= n>1)
    if(diag == TRUE){
        tmp <- cbind(seq(1, n), seq(0, n - 1))
    }
    if(diag == FALSE){
        tmp <- cbind(seq(1, n - 1), seq(0, n - 2))
    }
    tmp <- cbind(tmp, apply(tmp, 1, function(x) x[2] * n + x[1]))
    tmp <- cbind(tmp, apply(tmp, 1, function(x) x[1] * n))
    if(diag == TRUE){
        idx <- unlist(apply(tmp, 1, function(x) seq(x[3], x[4])))
    }
    if(diag == FALSE){
        idx <- unlist(apply(tmp, 1, function(x) seq(x[3] + 1, x[4])))
    }
    return(idx)
}
