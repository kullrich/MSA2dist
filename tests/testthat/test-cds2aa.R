data(hiv)

test_that("cds2aa() outputs AAStringSet", {
    expect_true(class(cds2aa(hiv)) == "AAStringSet")
    expect_true(as.character(cds2aa(Biostrings::DNAStringSet("ATGCAACATTGC")))
    == "MQHC")
    expect_true(length(cds2aa(Biostrings::DNAStringSet("ATGCAACATTGC"),
    frame = 2)) == 0)
    expect_true(as.character(cds2aa(Biostrings::DNAStringSet("ATGCAACATTGC"),
    frame = 2, shorten = TRUE)) == "CNI")
    cds1 <- Biostrings::DNAString("ATGCAACATTGC")
    cds2 <- Biostrings::DNAString("ATG---CATTGC")
    cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
        Biostrings::DNAStringSet(cds2))
    expect_true(as.character(cds2aa(cds1.cds2.aln,framelist=c(2,3),
    shorten=TRUE))[2] == "XXL")
    expect_true(as.character(cds2aa(Biostrings::DNAStringSet("GTG"),
    genetic.code=Biostrings::getGeneticCode("2"))) == "M")
})
