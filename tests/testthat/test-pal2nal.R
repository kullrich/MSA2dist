data(hiv)

test_that("pal2nal() outputs DNAStringSet", {
    cds1 <- Biostrings::DNAString("ATGCAACATTGC")
    cds2 <- Biostrings::DNAString("ATGCATTGC")
    cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
        Biostrings::DNAStringSet(cds2))
    names(cds1.cds2.aln) <- c("cds1", "cds2")
    aa1 <- Biostrings::AAString("MQHC")
    aa2 <- Biostrings::AAString("M-HC")
    aa1.aa2.aln <- c(Biostrings::AAStringSet(aa1),
        Biostrings::AAStringSet(aa2))
    names(aa1.aa2.aln) <- c("cds1", "cds2")
    expect_true(class(pal2nal(aa1.aa2.aln, cds1.cds2.aln)) == "DNAStringSet")
})
