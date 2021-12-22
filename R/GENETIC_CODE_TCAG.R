#' @title GENETIC_CODE_TCAG
#' @name GENETIC_CODE_TCAG
#' @description \code{GENETIC_CODE} from \code{Biostrings} extended by codon
#' number and number of syn sites.
#' @param codon \code{codon} [mandatory]
#' @return An object of class \code{numeric}
#' @importFrom Biostrings GENETIC_CODE
#' @importFrom stats setNames
#' @seealso \code{\link[Biostrings]{GENETIC_CODE}}
#' @examples
#' GENETIC_CODE_TCAG
#' @export GENETIC_CODE_TCAG
#' @author Kristian K Ullrich

codon2number <- function(codon){
    baseNumber <- setNames(c(0, 1, 2, 3), c("T", "C", "A", "G"))
    left <- baseNumber[str_sub(codon, 1, 1)]
    mid <- baseNumber[str_sub(codon, 2, 2)]
    right <- baseNumber[str_sub(codon, 3, 3)]
    codonnumber <- (mid * 16) + (left * 4) + (right) + 1
    names(codonnumber) <- codon
    return(codonnumber)
}

TMP_GENETIC_CODE <- Biostrings::GENETIC_CODE
TMP_CODON_NUMBER <- unlist(lapply(names(TMP_GENETIC_CODE), codon2number))
TMP_GENETIC_CODE <- TMP_GENETIC_CODE[order(TMP_CODON_NUMBER)]
TMP_CODON_NUMBER <- TMP_CODON_NUMBER[order(TMP_CODON_NUMBER)]
TMP_SYN_SITES <- setNames(c(1,1,2,2,
                            3,3,4,4,
                            2,2,2,0,
                            3,3,3,3,
                            3,3,3,3,
                            3,3,3,3,
                            3,3,3,3,
                            3,3,3,3,
                            1,1,2,1,
                            1,1,1,1,
                            1,1,1,1,
                            1,1,1,1,
                            1,1,1,0,
                            3,3,4,4,
                            1,1,2,2,
                            3,3,3,3), names(TMP_CODON_NUMBER))

GENETIC_CODE_TCAG <- data.frame(CODON = names(TMP_GENETIC_CODE),
                                GENETIC_CODE = TMP_GENETIC_CODE,
                                CODON_NUMBER = TMP_CODON_NUMBER,
                                SYN_SITES = TMP_SYN_SITES)
