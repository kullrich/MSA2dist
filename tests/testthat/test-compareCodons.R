test_that("compareCodons()", {
    for(i in names(Biostrings::GENETIC_CODE)){
    for(j in names(Biostrings::GENETIC_CODE)){
    compareCodons(i,j)}}
    expect_true(compareCodons(i,j)[1] == 0)
})
