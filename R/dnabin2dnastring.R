#' @title dnabin2dnastring
#' @name dnabin2dnastring
#' @description This function converts an \code{ape} \code{DNAbin} into a
#' \code{DNAStringSet}.
#' @param dnabin \code{ape} \code{DNAbin} [mandatory]
#' @return An object of class \code{DNAStringSet}
#' @importFrom methods is
#' @importFrom methods slot
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAString
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readAAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings width
#' @importFrom Biostrings subseq
#' @importFrom seqinr as.alignment
#' @importFrom ape as.DNAbin.alignment
#' @seealso \code{\link[seqinr]{as.alignment}}
#' \code{\link[ape]{as.DNAbin.alignment}}
#' \code{\link[Biostrings]{DNAStringSet}}
#' @examples
#' data(woodmouse, package="ape")
#' ## convert into DNAStringSet
#' #dnabin2dnastring(woodmouse)
#' woodmouse |> dnabin2dnastring()
#' @export dnabin2dnastring
#' @author Kristian K Ullrich

dnabin2dnastring <- function(dnabin){
    stopifnot("Error: input needs to be a DNAbin"=
        methods::is(dnabin, "DNAbin"))
    dna <- setNames(Biostrings::DNAStringSet(ape::as.alignment(dnabin)$seq),
        labels(dnabin))
    return(dna)
}
