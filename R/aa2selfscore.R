#' @title aa2selfscore
#' @name aa2selfscore
#' @description This function return the selfscore from an \code{AAStringSet}.
#' @param aa \code{AAStringSet} [mandatory]
#' @param scorematrix score matrix to use [default: BLOSUM62]
#' @return \code{data.frame}
#' @importFrom methods is slot
#' @importFrom Biostrings AAString AAStringSet
#' @importFrom stringr word
#' @importFrom utils data
#' @seealso \code{\link[Biostrings]{XStringSet-class}},
#' \code{\link[Biostrings]{substitution.matrices}}
#' @examples
#' data(woodmouse, package="ape")
#' #cds2aa(dnabin2dnastring(woodmouse), shorten=TRUE,
#' #genetic.code=Biostrings::getGeneticCode("2"))
#' woodmouse |> dnabin2dnastring() |> cds2aa(shorten=TRUE,
#' genetic.code=Biostrings::getGeneticCode("2")) |> aa2selfscore()
#' @export aa2selfscore
#' @author Kristian K Ullrich

aa2selfscore <- function(aa, scorematrix="BLOSUM62"){
    stopifnot("Error: input needs to be a AAStringSet"=
        methods::is(aa, "AAStringSet"))
    stopifnot("Error: score needs to be BLOSUM45, BLOSUM50, BLOSUM62,
        BLOSUM80, BLOSUM100, PAM30, PAM40, PAM70,
        PAM120 or PAM250"= scorematrix %in%
        c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30",
        "PAM40", "PAM70", "PAM120", "PAM250"))
    if(!is.null(names(aa))){
        names(aa) <- stringr::word(names(aa), 1)
    }
    BLOSUM45<-NULL
    BLOSUM50<-NULL
    BLOSUM62<-NULL
    BLOSUM80<-NULL
    BLOSUM100<-NULL
    PAM30<-NULL
    PAM40<-NULL
    PAM70<-NULL
    PAM120<-NULL
    PAM250<-NULL
    if(scorematrix == "BLOSUM45"){utils::data(BLOSUM45, package="Biostrings",
        envir=environment())
        score.matrix <- BLOSUM45}
    if(scorematrix == "BLOSUM50"){utils::data(BLOSUM50, package="Biostrings",
        envir=environment())
        score.matrix <- BLOSUM50}
    if(scorematrix == "BLOSUM62"){utils::data(BLOSUM62, package="Biostrings",
        envir=environment())
        score.matrix <- BLOSUM62}
    if(scorematrix == "BLOSUM80"){utils::data(BLOSUM80, package="Biostrings",
        envir=environment())
        score.matrix <- BLOSUM80}
    if(scorematrix == "BLOSUM100"){utils::data(BLOSUM100, package="Biostrings",
        envir=environment())
        score.matrix <- BLOSUM100}
    if(scorematrix == "PAM30"){utils::data(PAM30, package="Biostrings",
        envir=environment())
        score.matrix <- PAM30}
    if(scorematrix == "PAM40"){utils::data(PAM40, package="Biostrings",
        envir=environment())
        score.matrix <- PAM40}
    if(scorematrix == "PAM70"){utils::data(PAM70, package="Biostrings",
        envir=environment())
        score.matrix <- PAM70}
    if(scorematrix == "PAM120"){utils::data(PAM120, package="Biostrings",
        envir=environment())
        score.matrix <- PAM120}
    if(scorematrix == "PAM250"){utils::data(PAM250, package="Biostrings",
        envir=environment())
        score.matrix <- PAM250}
    aa_selfscore_list <- lapply(aa, function(x) {
        sum(unlist(lapply(strsplit(as.character(x), ""),
        function(y){diag(score.matrix)[y]})))})
    aa_selfscore <- data.frame(name=names(aa_selfscore_list),
        score=unlist(aa_selfscore_list))
    return(aa_selfscore)
}
