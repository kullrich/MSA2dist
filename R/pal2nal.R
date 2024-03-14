#' @title pal2nal
#' @name pal2nal
#' @description This function takes an \code{AAStringSet} alignment and
#' its corresponding coding sequences \code{DNAStringSet} and converts
#' the protein alignment into a codon alignment.
#' @param pal \code{AAStringSet} [mandatory]
#' @param nal \code{DNAStringSet} [mandatory]
#' @return codon alignment as \code{DNAStringSet}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom methods is slot
#' @references Pagès, H et al. (2014) Biostrings: Efficient manipulation of
#' biological strings. \emph{R package version}, \bold{2(0)}.
#' @seealso \code{\link[Biostrings]{pairwiseAlignment}}
#' @examples
#' ## define two cds sequences
#' cds <- Biostrings::DNAStringSet(c("ATGCAACATTGC", "ATGCATTGC"))
#' names(cds) <- c("cds1", "cds2")
#' ## get protein alignment
#' aa <- MSA2dist::cds2aa(cds)
#' msa <- makePostalignedSeqs(Biostrings::pairwiseAlignment(aa[1], aa[2]))[[1L]]
#' names(msa) <- names(aa)
#' ## get codon alignment
#' nal <- MSA2dist::pal2nal(pal=msa, nal=cds)
#' nal
#' @export pal2nal
#' @author Kristian K Ullrich

pal2nal <- function(pal,
    nal,
    remove.gaps=FALSE){
    nal <- nal[names(pal)]
    pal.gap.pos <- lapply(pal, function(x) {
        gregexpr("-+", x)})
    pal.gap.len <- lapply(pal.gap.pos, function(x) {
        attr(x[[1]], "match.length")})
    pal.gap.pos <- lapply(pal.gap.pos, function(x) {
        ifelse(x[[1]]!=-1, x[[1]], -1)})
    pal.gap.len <- lapply(pal.gap.len, function(x) {
        if (x[1] != -1) return(x) else return(0)})
    pal.gap.len.cumsum <- lapply(pal.gap.len, cumsum)
    nal_out <- Biostrings::DNAStringSet()
    for(i in seq(from=1, to=length(pal))){
        if(pal.gap.pos[[i]][1] == -1){
            nal_out <- c(nal_out, nal[i])
        } else {
            nal_i <- as.character(nal[[i]])
            n_i_codons <- nchar(nal_i)/3
            n_i <- ""
            n_i_codons_added <- 0
            for(j in seq(from=1, to=length(pal.gap.pos[[i]]))){
                gap_pos <- pal.gap.pos[[i]][j]
                gap_len <- pal.gap.len[[i]][j]
                gap_len_cumsum <- pal.gap.len.cumsum[[i]][j]
                gap <- paste0(rep("---", gap_len), collapse="")
                if(gap_pos==1){
                    n_i <- paste0(n_i, gap)
                } else {
                    n_i_codons_to_add <- gap_pos-n_i_codons_added
                    n_i<- paste0(n_i, substr(nal_i,
                        (n_i_codons_added*3)+1,
                        (n_i_codons_added+n_i_codons_to_add)*3))
                    n_i <- paste0(n_i, gap)
                    n_i_codons_added <- n_i_codons_added+n_i_codons_to_add
                }
            }
            if(n_i_codons_added!=n_i_codons){
                n_i<- paste0(n_i, substr(nal_i,
                    (n_i_codons_added*3)+1,
                    n_i_codons*3))
            }
            nal_out <- c(nal_out,
                setNames(Biostrings::DNAStringSet(n_i),
                names(nal)[i]))
        }
    }
    if(remove.gaps){
        nal_out <- Biostrings::DNAStringSet(apply(as.matrix(nal_out)[,
            apply(as.matrix(nal_out), 2, function(x) !any(x=="-"))], 1,
            function(x) paste(x, collapse="")))
        names(nal_out) <- names(nal)
    }
    return(nal_out)
}
