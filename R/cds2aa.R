#' @title cds2aa
#' @name cds2aa
#' @description This function translates a \code{DNAStringSet} into an
#' \code{AAStringSet}.
#' @param cds \code{DNAStringSet} [mandatory]
#' @param shorten shorten all sequences to multiple of three [default: FALSE]
#' @param frame  indicates the first base of a the first codon [default: 1]
#' @param framelist  supply vector of frames for each entry [default: NULL]
#' @param genetic.code The genetic code to use for the translation of codons
#' into Amino Acid letters [default: NULL]
#' @return \code{AAStringSet}
#' @importFrom methods is slot
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq translate
#' getGeneticCode GENETIC_CODE GENETIC_CODE_TABLE
#' @importFrom stringr word
#' @seealso \code{\link[Biostrings]{XStringSet-class}},
#' \code{\link[seqinr]{translate}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' #cds2aa(cds1.cds2.aln)
#' cds1.cds2.aln |> cds2aa()
#' ## alternative genetic code
#' data(woodmouse, package="ape")
#' #cds2aa(dnabin2dnastring(woodmouse), shorten=TRUE)
#' woodmouse |> dnabin2dnastring() |> cds2aa(shorten=TRUE)
#' #cds2aa(dnabin2dnastring(woodmouse), shorten=TRUE,
#' woodmouse |> dnabin2dnastring() |> cds2aa(shorten=TRUE,
#' genetic.code=Biostrings::getGeneticCode("2"))
#' @export cds2aa
#' @author Kristian K Ullrich

cds2aa <- function(cds, shorten=FALSE, frame=1, framelist=NULL,
    genetic.code=NULL){
    stopifnot("Error: input needs to be a DNAStringSet"=
        methods::is(cds, "DNAStringSet"))
    stopifnot("Error: frame needs to be 1 or 2 or 3"= frame %in% c(1, 2, 3))
    if(!is.null(framelist)){
        stopifnot("Error: framelist needs to be of equal length as cds"=
            length(framelist) == length(cds))
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
    cds <- Biostrings::DNAStringSet(gsub("-", "N", cds))
    cds <- Biostrings::DNAStringSet(gsub("X", "N", cds))
    if(!is.null(genetic.code)){
        aa <- Biostrings::translate(cds, genetic.code=genetic.code,
                                    if.fuzzy.codon="X")
        return(aa)
    }
    aa <- Biostrings::translate(cds, if.fuzzy.codon="X")
    return(aa)
}
