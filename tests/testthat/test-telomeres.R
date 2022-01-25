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
