#' @title iupacMatrix
#' @name iupacMatrix
#' @description This function creates a \code{iupacMatrix} object to be used
#' with the \code{rcpp_distSTRING} function. By default, the \code{iupac matrix}
#' is defined as literal distance obtained from \code{Chang et al. 2017}.
#' (see \url{https://link.springer.com/article/10.1007/s00335-017-9704-9})
#' @return \code{score matrix}
#' @references Chang,P. L.,Kopania,E.,Keeble,S.,Sarver,B. A.,Larson,
#' E.,Orth,A.,... & Dean,M. D. (2017). Whole exome sequencing of
#' wild-derived inbred strains of mice improves power to link phenotype and
#' genotype. \emph{Mammalian genome},\bold{28(9-10)},416-425.
#' @seealso \link[MSA2dist]{dnastring2dist},\link[ape]{dist.dna}
#' @examples
#' iupacMatrix()
#' @export iupacMatrix
#' @author Kristian K Ullrich
## A C G T
## R Y S W K M
## B D H V
## . - N X

iupacMatrix<-function(){
    distances<-c(
    # A
    0,1,1,1,0.5,1,1,0.5,1,0.5,-1,-1,-1,-1,-1,-1,-1,-1,
    # C
    1,0,1,1,1,0.5,0.5,1,1,0.5,-1,-1,-1,-1,-1,-1,-1,-1,
    # G
    1,1,0,1,0.5,1,0.5,1,0.5,1,-1,-1,-1,-1,-1,-1,-1,-1,
    # T
    1,1,1,0,1,0.5,1,0.5,0.5,1,-1,-1,-1,-1,-1,-1,-1,-1,
    # R
    0.5,1,0.5,1,0,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,
    # Y
    1,0.5,1,0.5,1,0,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,
    # S
    1,0.5,0.5,1,1,1,0,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,
    # W
    0.5,1,1,0.5,1,1,1,0,1,1,-1,-1,-1,-1,-1,-1,-1,-1,
    # K
    1,1,0.5,0.5,1,1,1,1,0,1,-1,-1,-1,-1,-1,-1,-1,-1,
    # M
    0.5,0.5,1,1,1,1,1,1,1,0,-1,-1,-1,-1,-1,-1,-1,-1,
    # B
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    # D
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    # H
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    # V
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    # .
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    # -
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    # N
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    # X
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
    )
    iupacMatrix<-matrix(distances,ncol=18,nrow=18)
    colnames(iupacMatrix)<-rownames(iupacMatrix)<-c("A","C","G","T","R","Y",
    "S","W","K","M","B","D","H","V",".","-","N","X")
    return(iupacMatrix)
}
