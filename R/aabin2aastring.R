#' @title aabin2aastring
#' @name aabin2aastring
#' @description This function converts an \code{ape} \code{AAbin} into
#' \code{AAStringSet}.
#' @param aabin \code{ape} \code{AAbin} [mandatory]
#' @return An object of class \code{AAStringSet}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom seqinr as.alignment
#' @importFrom ape as.DNAbin.alignment
#' @seealso \code{\link[seqinr]{as.alignment}}
#' \code{\link[ape]{as.DNAbin.alignment}}
#' \code{\link[Biostrings]{AAStringSet}}
#' @examples
#' data(woodmouse, package="ape")
#' ## convert into AAStringSet
#' #aabin2aastring(ape::trans(woodmouse, 2))
#' ape::trans(woodmouse, 2) |> aabin2aastring()
#' @export aabin2aastring
#' @author Kristian K Ullrich

aabin2aastring <- function(aabin){
    stopifnot("Error: input needs to be an AAbin"=
        methods::is(aabin, "AAbin"))
    if(!is.null(dim(aabin))){
        aa <- setNames(Biostrings::AAStringSet(unlist(apply(
            as.character(aabin), 1, function(x) paste0(x, collapse = "")))),
            labels(aabin))
    }
    if(is.null(dim(aabin))){
        aa <- setNames(Biostrings::AAStringSet(unlist(lapply(
            as.character(aabin), function(x) paste0(x, collapse = "")))),
            labels(aabin))
    }
    return(aa)
}
