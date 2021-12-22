#' @title subString
#' @name subString
#' @description This function gets a subsequence from a \code{DNAString},
#' \code{RNAString}, \code{AAString}, \code{BString}, \code{DNAStringSet},
#' \code{RNAStringSet}, \code{AAStringSet}, \code{BStringSet}
#' object from the \code{Biostrings} package.
#' @return subsequence of an \code{Biostrings} object
#' @importFrom Biostrings DNAString RNAString AAString BString DNAStringSet
#' RNAStringSet AAStringSet BStringSet
#' @param x \code{DNAStringSet}, \code{RNAString}, \code{AAString},
#' \code{BString}, \code{DNAStringSet}, \code{RNAStringSet},
#' \code{AAStringSet}, \code{BStringSet} [mandatory]
#' @param s start vector [mandatory]
#' @param e end vector [mandatory]
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'     Biostrings::DNAStringSet(cds2))
#' subString(cds1.cds2.aln, c(1,7), c(3,12))
#' @seealso \link[XVector]{subseq}
#' @export subString
#' @author Kristian K Ullrich

subString<-function(x, s, e){
    x.class<-class(x)[1]
    se.matrix<-cbind(s, e)
    x.len <- length(x)
    if(x.class=="DNAString"){
        myDNAString<-x
        newDNAString<-Biostrings::DNAString(paste(
            vapply(
                apply(se.matrix, 1, function(x) {
                    Biostrings::subseq(myDNAString, x[1], x[2])
                }), function(x) {
                paste0(x)
            }, rep("", 1)), sep="", collapse=""))
        names(newDNAString)<-names(x)
        return(newDNAString)
    }
    if(x.class=="RNAString"){
        myRNAString<-x
        newRNAString<-Biostrings::RNAString(paste(
            vapply(
                apply(se.matrix, 1, function(x) {
                    Biostrings::subseq(myRNAString, x[1], x[2])
                }), function(x) {
                paste0(x)
            }, rep("", 1)), sep="", collapse=""))
        names(newRNAString)<-names(x)
        return(newRNAString)
    }
    if(x.class=="AAString"){
        myAAString<-x
        newAAString<-Biostrings::AAString(paste(
            vapply(
                apply(se.matrix, 1, function(x) {
                    Biostrings::subseq(myAAString, x[1], x[2])
                }), function(x) {
                paste0(x)
            }, rep("", 1)), sep="", collapse=""))
        names(newAAString)<-names(x)
        return(newAAString)
    }
    if(x.class=="BString"){
        myBString<-x
        newBString<-Biostrings::BString(paste(
            vapply(
                apply(se.matrix, 1, function(x) {
                    Biostrings::subseq(myBString, x[1], x[2])
                }), function(x) {
                paste0(x)
            }, rep("", 1)), sep="", collapse=""))
        names(newBString)<-names(x)
        return(newBString)
    }
    if(x.class=="DNAStringSet"){
        myDNAStringSet<-x
        if(length(myDNAStringSet)>1){
            newDNAStringSet<-Biostrings::DNAStringSet(
                apply(
                    vapply(
                        apply(se.matrix, 1, function(x) {
                            Biostrings::subseq(myDNAStringSet, x[1], x[2])
                        }), function(x) {
                        paste0(x)
                    }, rep("", length(myDNAStringSet))), 1, function(x) {
                    paste(x, sep="", collapse="")
                }))
            names(newDNAStringSet)<-names(x)
            return(newDNAStringSet)
        }
        if(length(myDNAStringSet)==1){
            newDNAStringSet<-Biostrings::DNAStringSet(paste(
                vapply(
                    apply(se.matrix, 1, function(x) {
                        Biostrings::subseq(myDNAStringSet, x[1], x[2])
                    }), function(x) {
                    paste0(x)
                }, rep("", length(myDNAStringSet))), sep="", collapse=""))
            names(newDNAStringSet) <- names(x)
            return(newDNAStringSet)
        }
    }
    if(x.class=="RNAStringSet"){
        myRNAStringSet<-x
        if(length(myRNAStringSet)>1){
            newRNAStringSet<-Biostrings::RNAStringSet(
                apply(
                    vapply(
                        apply(se.matrix, 1, function(x) {
                            Biostrings::subseq(myRNAStringSet, x[1], x[2])
                        }), function(x) {paste0(x)
                    }, rep("", length(myRNAStringSet))), 1, function(x) {
                    paste(x, sep="", collapse="")
                }))
            names(newRNAStringSet)<-names(x)
            return(newRNAStringSet)
        }
        if(length(myRNAStringSet)==1){
            newRNAStringSet<-Biostrings::RNAStringSet(paste(
                vapply(
                    apply(se.matrix, 1, function(x) {
                        Biostrings::subseq(myRNAStringSet, x[1], x[2])
                    }), function(x) {
                    paste0(x)
                }, rep("", length(myRNAStringSet))), sep="", collapse=""))
            names(newRNAStringSet)<-names(x)
            return(newRNAStringSet)
        }
    }
    if(x.class=="AAStringSet"){
        myAAStringSet<-x
        if(length(myAAStringSet)>1){
            newAAStringSet<-Biostrings::AAStringSet(
                apply(
                    vapply(
                        apply(se.matrix, 1, function(x) {
                            Biostrings::subseq(myAAStringSet, x[1], x[2])
                        }), function(x) {
                        paste0(x)
                    }, rep("", length(myAAStringSet))), 1, function(x) {
                    paste(x, sep="", collapse="")
                }))
            names(newAAStringSet)<-names(x)
            return(newAAStringSet)
        }
        if(length(myAAStringSet)==1){
            newAAStringSet<-Biostrings::AAStringSet(paste(
                vapply(
                    apply(se.matrix, 1, function(x) {
                        Biostrings::subseq(myAAStringSet, x[1], x[2])
                    }), function(x) {
                    paste0(x)
                }, rep("", length(myAAStringSet))), sep="", collapse=""))
            names(newAAStringSet)<-names(x)
            return(newAAStringSet)
        }
    }
    if(x.class=="BStringSet"){
        myBStringSet<-x
        if(length(myBStringSet)>1){
            newBStringSet<-Biostrings::BStringSet(
                apply(
                    vapply(
                        apply(se.matrix, 1, function(x) {
                            Biostrings::subseq(myBStringSet, x[1], x[2])
                        }), function(x) {
                        paste0(x)
                    }, rep("", length(myBStringSet))), 1, function(x) {
                    paste(x, sep="", collapse="")
                }))
            names(newBStringSet)<-names(x)
            return(newBStringSet)
        }
        if(length(myBStringSet)==1){
            newBStringSet<-Biostrings::BStringSet(paste(
                vapply(
                    apply(se.matrix, 1, function(x) {
                        Biostrings::subseq(myBStringSet, x[1], x[2])
                        }), function(x) {
                        paste0(x)
                }, rep("", length(myBStringSet))), sep="", collapse=""))
            names(newBStringSet)<-names(x)
            return(newBStringSet)
        }
    }
}
