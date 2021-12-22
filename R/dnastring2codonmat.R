#' @title dnastring2codonmat
#' @name dnastring2codonmat
#' @description This function converts a \code{DNAStringSet} into a
#' \code{codon matrix}.
#' @param cds \code{DNAStringSet} [mandatory]
#' @param shorten shorten all sequences to multiple of three [default: FALSE]
#' @param frame  indicates the first base of a the first codon [default: 1]
#' @param framelist  supply vector of frames for each entry [default: NULL]
#' @return An object of class \code{alignment} which is a list with the
#' following components:\cr
#' \code{nb} the number of aligned sequences\cr
#' \code{nam} a vector of strings containing the names of the aligned
#' sequences\cr
#' \code{seq} a vector of strings containing the aligned sequences\cr
#' \code{com} a vector of strings containing the commentaries for each sequence
#' or \code{NA} if there are no comments
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom stringr word
#' @seealso \code{\link[seqinr]{as.alignment}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' ## convert into alignment
#' #dnastring2codonmat(cds1.cds2.aln)
#' cds1.cds2.aln |> dnastring2codonmat()
#' ## use frame 2 and shorten to circumvent multiple of three error
#' cds1 <- Biostrings::DNAString("-ATGCAACATTGC-")
#' cds2 <- Biostrings::DNAString("-ATG---CATTGC-")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' cds1.cds2.aln |> dnastring2codonmat(frame=2, shorten=TRUE)
#' @export dnastring2codonmat
#' @author Kristian K Ullrich

dnastring2codonmat <- function(cds, shorten = FALSE, frame = 1, framelist=NULL){
    stopifnot("Error: input needs to be a DNAStringSet"=
        methods::is(cds, "DNAStringSet"))
    stopifnot("Error: input needs to be a Alignment of equal width"=
        length(unique(Biostrings::width(cds)))==1)
    stopifnot("Error: frame needs to be 1 or 2 or 3"= frame %in% c(1, 2, 3))
    if(!is.null(framelist)){
        stopifnot("Error: framelist needs to be of equal length as cds"=
            length(framelist) != length(cds))
    }
    if(!is.null(names(cds))){
        names(cds) <- stringr::word(names(cds), 1)
    }
    if(is.null(framelist)){
        cds <- Biostrings::subseq(cds, frame, Biostrings::width(cds))
    }
    if(!is.null(framelist)){
        cds <- Biostrings::subseq(cds, framelist, Biostrings::width(cds))
    }
    if(shorten){
        cds <- Biostrings::subseq(cds, 1,
            Biostrings::width(cds) - Biostrings::width(cds) %% 3)
    }
    cds_not_multiple_of_three.idx <- which(Biostrings::width(cds) %% 3 != 0)
    if(length(cds_not_multiple_of_three.idx) > 0){
        cds_not_multiple_of_three <- cds[cds_not_multiple_of_three.idx]
        cds <- cds[-cds_not_multiple_of_three.idx]
    }
    # remove IUPAC chars
    cds <- Biostrings::DNAStringSet(gsub("R", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("Y", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("S", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("W", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("K", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("M", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("B", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("D", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("H", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("V", "N", cds))
    cds.codonmat <- apply(cbind(Biostrings::width(cds), as.character(cds)), 1,
        function(x) {stringi::stri_sub(x[2], seq(1, as.numeric(x[1]), by = 3),
        length = 3)})
    return(cds.codonmat)
}
