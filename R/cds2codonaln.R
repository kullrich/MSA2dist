#' @title cds2codonaln
#' @name cds2codonaln
#' @description This function takes two single sequence \code{DNAString}'s or
#' two single sequence \code{DNAStringSet}'s, converts them into aa, calculates
#' a global alignment and converts this alignment back into a codon alignment.
#' @param cds1 single sequence \code{DNAStringSet} or \code{DNAString}
#' [mandatory]
#' @param cds2 single sequence \code{DNAStringSet} or \code{DNAString}
#' [mandatory]
#' @param type type of alignment (see
#' \code{\link[Biostrings]{pairwiseAlignment}}) [default: global]
#' @param substitutionMatrix substitution matrix representing the fixed
#' substitution scores for an alignment (see
#' \code{\link[Biostrings]{pairwiseAlignment}}) [default: BLOSUM62]
#' @param gapOpening the cost for opening a gap in the alignment (see
#' \code{\link[Biostrings]{pairwiseAlignment}}) [default: 10]
#' @param gapExtension the incremental cost incurred along the length of the
#' gap in the alignment (see \code{\link[Biostrings]{pairwiseAlignment}})
#' [default: 0.5]
#' @param remove.gaps specify if gaps in the codon alignment should be removed
#' [default: FALSE]
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
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATGCATTGC")
#' cds2codonaln(cds1, cds2)
#' @export cds2codonaln
#' @author Kristian K Ullrich

cds2codonaln <- function(cds1, cds2,
    type="global",
    substitutionMatrix="BLOSUM62",
    gapOpening=10,
    gapExtension=0.5,
    remove.gaps=FALSE){
    stopifnot("Error: cds1 needs to be either DNAString or DNAStringSet"=
        {methods::is(cds1, "DNAString") || methods::is(cds1, "DNAStringSet")})
    stopifnot("Error: cds2 needs to be either DNAString or DNAStringSet"=
        {methods::is(cds2, "DNAString") || methods::is(cds2, "DNAStringSet")})
    if(methods::is(cds1, "DNAString")){
        x.aa <- MSA2dist::cds2aa(Biostrings::DNAStringSet(cds1))[[1]]
        x.name <- "cds1"
    }
    if(methods::is(cds2, "DNAString")){
        y.aa <- MSA2dist::cds2aa(Biostrings::DNAStringSet(cds2))[[1]]
        y.name <- "cds2"
    }
    if(methods::is(cds1, "DNAStringSet")){
        stopifnot("Error: cds1 needs to only contain one sequence"=
            length(cds1) == 1)
        x.aa <- MSA2dist::cds2aa(cds1)[[1]]
        x.name <- names(cds1)
        cds1 <- cds1[[1]]
    }
    if(methods::is(cds2, "DNAStringSet")){
        stopifnot("Error: cds2 needs to only contain one sequence"=
            length(cds2) == 1)
        y.aa <- MSA2dist::cds2aa(cds2)[[1]]
        y.name <- names(cds2)
        cds2 <- cds2[[1]]
    }
    xy.aln <- makePostalignedSeqs(Biostrings::pairwiseAlignment(x.aa, y.aa,
        type=type, substitutionMatrix=substitutionMatrix, gapOpening=gapOpening,
        gapExtension=gapExtension))[[1L]]
    xy.aln.pattern.gap.pos <- gregexpr("-+", xy.aln[1])
    xy.aln.pattern.gap.len <- attr(xy.aln.pattern.gap.pos[[1]], "match.length")
    if(xy.aln.pattern.gap.pos[[1]][1]!=-1){
        xy.aln.pattern.gap.pos <- xy.aln.pattern.gap.pos[[1]]
        xy.aln.pattern.gap.len <- xy.aln.pattern.gap.len
    }
    if(xy.aln.pattern.gap.pos[[1]][1]==-1){
        xy.aln.pattern.gap.pos <- -1
        xy.aln.pattern.gap.len <- 0
    }
    xy.aln.subject.gap.pos <- gregexpr("-+", xy.aln[2])
    xy.aln.subject.gap.len <- attr(xy.aln.subject.gap.pos[[1]], "match.length")
    if(xy.aln.subject.gap.pos[[1]][1]!=-1){
        xy.aln.subject.gap.pos <- xy.aln.subject.gap.pos[[1]]
        xy.aln.subject.gap.len <- xy.aln.subject.gap.len
    }
    if(xy.aln.subject.gap.pos[[1]][1]==-1){
        xy.aln.subject.gap.pos <- -1
        xy.aln.subject.gap.len <- 0
    }
    x.cds <- Biostrings::subseq(cds1, (xy.aln@ranges@start[1]*3)-2,
        (xy.aln@ranges@width[1]-sum(xy.aln.pattern.gap.len))*3)
    y.cds <- Biostrings::subseq(cds2, (xy.aln@ranges@start[2]*3)-2,
        (xy.aln@ranges@width[2]-sum(xy.aln.subject.gap.len))*3)
    tmp.x <- Biostrings::DNAString()
    tmp.x.cur.start <- 0
    tmp.x.cur.end <- 0
    for(i in seq(from=1, to=(length(xy.aln.pattern.gap.pos)))){
        cur.gap.pos <- xy.aln.pattern.gap.pos[i]
        cur.gap.len <- xy.aln.pattern.gap.len[i]
        if(i==1){
            if(cur.gap.pos==1){
                #add gaps
                gaps <- Biostrings::DNAString(paste(rep("---", cur.gap.len),
                    collapse=""))
                tmp.x <- c(tmp.x, gaps)
            }
            if(cur.gap.pos!=1){
                #add seq to first gap position
                tmp.x.cur.start <- 1
                tmp.x.cur.end <- (cur.gap.pos-1)*3
                tmp.x <- c(tmp.x, Biostrings::subseq(x.cds, tmp.x.cur.start,
                    tmp.x.cur.end))
                #add gaps
                gaps <- Biostrings::DNAString(paste(rep("---", cur.gap.len),
                    collapse=""))
                tmp.x <- c(tmp.x, gaps)
            }
        }
        if(i!=1){
            tmp.x.cur.start <- tmp.x.cur.end+1
            tmp.x.cur.end <- tmp.x.cur.end+
                (cur.gap.pos-xy.aln.pattern.gap.pos[i-1]-
                xy.aln.pattern.gap.len[i-1])*3
            tmp.x<-c(tmp.x, Biostrings::subseq(x.cds, tmp.x.cur.start,
                tmp.x.cur.end))
            #add gaps
            gaps <- Biostrings::DNAString(paste(rep("---",cur.gap.len),
                collapse=""))
            tmp.x <- c(tmp.x, gaps)
        }
    }
    if(tmp.x.cur.end!=length(x.cds)){
        tmp.x.cur.start <- tmp.x.cur.end+1
        tmp.x.cur.end <- length(x.cds)
        tmp.x <- c(tmp.x, Biostrings::subseq(x.cds, tmp.x.cur.start,
            tmp.x.cur.end))
    }
    tmp.y <- Biostrings::DNAString()
    tmp.y.cur.start <- 0
    tmp.y.cur.end <- 0
    for(i in seq(from=1, to=(length(xy.aln.subject.gap.pos)))){
        cur.gap.pos<-xy.aln.subject.gap.pos[i]
        cur.gap.len<-xy.aln.subject.gap.len[i]
        if(i == 1){
            if(cur.gap.pos==1){
                #add gaps
                gaps <- Biostrings::DNAString(paste(rep("---", cur.gap.len),
                    collapse=""))
                tmp.y <- c(tmp.y, gaps)
            }
            if(cur.gap.pos!=1){
                #add seq to first gap position
                tmp.y.cur.start <- 1
                tmp.y.cur.end <- (cur.gap.pos-1)*3
                tmp.y <- c(tmp.y, Biostrings::subseq(y.cds, tmp.y.cur.start,
                    tmp.y.cur.end))
                #add gaps
                gaps <- Biostrings::DNAString(paste(rep("---", cur.gap.len),
                    collapse=""))
                tmp.y <- c(tmp.y, gaps)
            }
        }
        if(i!=1){
            tmp.y.cur.start <- tmp.y.cur.end+1
            tmp.y.cur.end <- tmp.y.cur.end+
                (cur.gap.pos-xy.aln.subject.gap.pos[i-1]-
                xy.aln.subject.gap.len[i-1])*3
            tmp.y <- c(tmp.y, Biostrings::subseq(y.cds, tmp.y.cur.start,
                tmp.y.cur.end))
            #add gaps
            gaps <- Biostrings::DNAString(paste(rep("---", cur.gap.len),
                collapse=""))
            tmp.y <- c(tmp.y, gaps)
        }
    }
    if(tmp.y.cur.end!=length(y.cds)){
        tmp.y.cur.start <- tmp.y.cur.end+1
        tmp.y.cur.end <- length(y.cds)
        tmp.y <- c(tmp.y, Biostrings::subseq(y.cds, tmp.y.cur.start,
            tmp.y.cur.end))
    }
    xy.cds.aln <- c(Biostrings::DNAStringSet(tmp.x),
        Biostrings::DNAStringSet(tmp.y))
    names(xy.cds.aln) <- c(x.name, y.name)
    if(remove.gaps){
        xy.cds.aln <- Biostrings::DNAStringSet(apply(as.matrix(xy.cds.aln)[,
            apply(as.matrix(xy.cds.aln), 2, function(x) !any(x=="-"))], 1,
            function(x) paste(x, collapse="")))
        names(xy.cds.aln) <- c(x.name, y.name)
    }
    return(xy.cds.aln)
}
