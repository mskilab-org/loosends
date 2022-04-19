library(loosends)

## unit test telomere detection
test_that(desc = "check C/G telomere detection runs", code = {
    suppressWarnings(
        expr = {
            test.sequences = c("CACCCTGACCCCAACCCT", "AGTTGGGGTTGGGGTTGGG", "AAAAAAAAAAA")
            gtel = find_telomeres(seq = test.sequences, gorc = "g")
            ctel = find_telomeres(seq = test.sequences, gorc = "c")
            expect_true(all(gtel == c(FALSE, TRUE, FALSE)))
            expect_true(all(ctel == c(TRUE, FALSE, FALSE)))
        }
    )
})

## unit test telomere detection
test_that(desc = "check C/G telomere detection runs", code = {
    suppressWarnings(
        expr = {
            test.sequences = c("TTAGGGTTAGGGTGAGGG",
                               "TTAGGGTTAGGGTTAGGGTT",
                               "CCCTAACCCTAACCCTAACC",
                               "CCGTAACCCTAACCATAA")
            dt = find_telomeres2(seq = test.sequences, gorc = "g")
            expect_true(all(dt[, grtr_canonical] == c(FALSE, TRUE, FALSE, FALSE)))
            expect_true(all(dt[, crtr_canonical] == c(FALSE, FALSE, TRUE, FALSE)))
            expect_true(all(dt[, grtr_noncanonical] == c(TRUE, FALSE, FALSE, FALSE)))
            expect_true(all(dt[, crtr_noncanonical] == c(FALSE, FALSE, FALSE, TRUE)))
        }
    )
})
