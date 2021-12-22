#' @title codon2numberAMBIG
#' @name codon2numberAMBIG
#' @description This function converts a \code{codon} into a \code{number},
#' but accept N and -.
#' @param codon [mandatory]
#' @return An object of class \code{numeric}
#' @importFrom Biostrings GENETIC_CODE
#' @importFrom stringr str_sub
#' @importFrom stats setNames
#' @seealso \code{\link[Biostrings]{GENETIC_CODE}}
#' @examples
#' #unlist(lapply(names(Biostrings::GENETIC_CODE), codon2numberAMBIG))
#' names(Biostrings::GENETIC_CODE) |> codon2numberAMBIG()
#' @export codon2numberAMBIG
#' @author Kristian K Ullrich

codon2numberAMBIG <- function(codon){
    baseNumber <- setNames(c(0, 1, 2, 3, 4, 5), c("T", "C", "A", "G", "N", "-"))
    left <- baseNumber[stringr::str_sub(codon, 1, 1)]
    mid <- baseNumber[stringr::str_sub(codon, 2, 2)]
    right <- baseNumber[stringr::str_sub(codon, 3, 3)]
    codonnumber <- (mid * 36) + (left * 6) + (right) + 1
    names(codonnumber) <- codon
    return(codonnumber)
}
