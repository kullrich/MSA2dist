#' @title codonmat2pnps
#' @name codonmat2pnps
#' @description This function calculates pn/ps according to \emph{Nei and
#' Gojobori (1986)}.
#' @param codonmat \code{codon matrix} of two columns to be
#' compared [mandatory]
#' @return An object of class \code{pnps} which is a list with the following
#' components:\cr
#' \code{seq1} sequence1 name\cr
#' \code{seq2} sequence2 name\cr
#' \code{Codons} sequence2 name\cr
#' \code{Compared} sequence2 name\cr
#' \code{Ambigiuous} sequence2 name\cr
#' \code{Indels} sequence2 name\cr
#' \code{Ns} sequence2 name\cr
#' \code{Sd} sequence2 name\cr
#' \code{Sn} sequence2 name\cr
#' \code{S} sequence2 name\cr
#' \code{N} sequence2 name\cr
#' \code{ps} sequence2 name\cr
#' \code{pn} sequence2 name\cr
#' \code{pnps} sequence2 name\cr
#' \code{ds} sequence2 name\cr
#' \code{dn} sequence2 name\cr
#' \code{dnds} sequence2 name\cr
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom seqinr kaks
#' @importFrom stats setNames
#' @seealso \code{\link[seqinr]{kaks}}
#' @references Nei and Gojobori. (1986) Simple methods for estimating the
#' numbers of synonymous and nonsynonymous nucleotide substitutions.
#' \emph{Mol. Biol. Evol.}, \bold{3(5)}, 418-426.
#' @references Ganeshan et al. (1997) Human immunodeficiency virus type 1
#' genetic evolution in children with different rates of development of disease.
#' \emph{J. Virology.} \bold{71(1)}, 663-677.
#' @references Yang et al. (2000) Codon-substitution models for heterogeneous
#' selection pressure at amino acid sites. \emph{Genetics.}
#' \bold{155(1)}, 431-449.
#' @examples
#' ## load example sequence data
#' data("hiv", package="MSA2dist")
#' #codonmat2pnps(dnastring2codonmat(hiv)[,c(1, 2)])
#' (hiv |> dnastring2codonmat())[,c(1, 2)] |> codonmat2pnps()
#' @export codonmat2pnps
#' @author Kristian K Ullrich

codonmat2pnps <- function(codonmat){
    stopifnot("Error: input needs to be a codonmat of 2 columns"=
        dim(codonmat)[2] == 2)
    seq1_name <- colnames(codonmat)[1]
    seq2_name <- colnames(codonmat)[2]
    count_codons <- dim(codonmat)[1]
    insertions_idx <- unique(unlist(apply(codonmat, 2,
        function(x) grep("-", x))))
    count_insertions <- length(insertions_idx)
    if(count_insertions > 0){
        codonmat <- codonmat[-insertions_idx, , drop = FALSE]
    }
    codonnumber <- apply(codonmat, 2, MSA2dist::codon2numberTCAG)
    Ns_idx <- unique(unlist(apply(codonnumber, 2, function(x) which(is.na(x)))))
    count_Ns <- length(Ns_idx)
    if(count_Ns > 0){
        codonmat <- codonmat[-Ns_idx, , drop = FALSE]
        codonnumber <- codonnumber[-Ns_idx, , drop = FALSE]
    }
    SA_Nei <- sum(MSA2dist::GENETIC_CODE_TCAG[codonmat[, 1], 4])
    SB_Nei <- sum(MSA2dist::GENETIC_CODE_TCAG[codonmat[, 2], 4])
    identical_codons_idx <- which(codonnumber[, 1] == codonnumber[, 2])
    identical_codons <- length(identical_codons_idx)
    if(identical_codons > 0){
        codonmat <- codonmat[-identical_codons_idx, , drop = FALSE]
        codonnumber <- codonnumber[-identical_codons_idx, , drop = FALSE]
    }
    ## At this point, codonA and codonB are "real" codons (no N's or -'s)
    ## but are not identical
    syn_codons <- 0
    nonsyn_codons <- 0
    if(nrow(codonmat) > 0){
        syn_nonsyn_codons <- apply(codonmat, 1,
            function(x) MSA2dist::compareCodons(x[1], x[2]))
        syn_codons <- syn_codons + sum(syn_nonsyn_codons[1, ])
        nonsyn_codons <- nonsyn_codons + sum(syn_nonsyn_codons[2, ])
    }
    count_ambiguous_codons <- count_insertions + count_Ns
    count_compared_codons <- count_codons - count_ambiguous_codons
    potential_syn <- ((SA_Nei / 3) + (SB_Nei / 3)) / 2
    potential_nonsyn <- (3 * count_compared_codons) - potential_syn
    ps <- syn_codons / potential_syn
    pn <- nonsyn_codons / potential_nonsyn
    ds <- ifelse(ps < 0.75 , ((-3 / 4) * log(1 - (4 * (ps / 3)))), NA)
    dn <- ifelse(pn < 0.75, ((-3 / 4) * log(1 - (4 * (pn / 3)))), NA)
    dnds <- NA
    if(dn != 0 & ds != 0 & !is.na(dn) & !is.na(ds)){ dnds <- dn/ds}
    pnps <- NA
    if(pn != 0 & ps != 0 & !is.na(pn) & !is.na(ps)){ pnps <- pn/ps}
    codonmat_out <- setNames(c(seq1_name, seq2_name, count_codons,
        count_compared_codons, count_ambiguous_codons, count_insertions,
        count_Ns, syn_codons, nonsyn_codons, potential_syn, potential_nonsyn,
        ps, pn, pnps, ds, dn, dnds),
        c("seq1", "seq2", "Codons", "Compared", "Ambigiuous", "Indels", "Ns",
        "Sd", "Sn", "S", "N", "ps", "pn", "pn/ps", "ds", "dn", "dn/ds"))
    attr(codonmat_out, "class") <- "pnps"
    return(codonmat_out)
}
