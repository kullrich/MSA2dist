#' @title codon2numberTCAG
#' @name codon2numberTCAG
#' @description This function converts a \code{codon} into a \code{number}.
#' @param codon [mandatory]
#' @return An object of class \code{numeric}
#' @importFrom Biostrings GENETIC_CODE
#' @importFrom stringr str_sub
#' @importFrom stats setNames
#' @seealso \code{\link[Biostrings]{GENETIC_CODE}}
#' @examples
#' #unlist(lapply(names(Biostrings::GENETIC_CODE), codon2numberTCAG))
#' names(Biostrings::GENETIC_CODE) |> codon2numberTCAG()
#' @export codon2numberTCAG
#' @author Kristian K Ullrich

codon2numberTCAG <- function(codon){
    baseNumber <- setNames(c(0, 1, 2, 3), c("T", "C", "A", "G"))
    left <- baseNumber[stringr::str_sub(codon, 1, 1)]
    mid <- baseNumber[stringr::str_sub(codon, 2, 2)]
    right <- baseNumber[stringr::str_sub(codon, 3, 3)]
    codonnumber <- (mid * 16) + (left * 4) + (right) + 1
    names(codonnumber) <- codon
    return(codonnumber)
}
