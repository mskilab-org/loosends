library(loosends)
library(testthat)

## TEST CASE 1
## files for first test case
loosereads.bam = system.file("tests", "testthat", "data", "loosereads_1", "loosereads.bam", package = "loosends")
aln.bam = system.file("tests", "testthat", "data", "loosereads_1", "aln.bam", package = "loosends")
qnames.txt = system.file("tests", "testthat", "data", "loosereads_1", "qnames.txt", package = "loosends")
windows.bed = system.file("tests", "testthat", "data", "loosereads_1", "windows.bed", package = "loosends")

## loci for first test case
this.le = parse.gr("20:60158837-")

## expected results for the first test case
qnames = readLines(con = qnames.txt)
windows = rtracklayer::import.bed(con = windows.bed)

test_that(desc = "loosereads_params", code = {
    suppressWarnings(
        expr = {
            this.params = grab_looseread_params(gr = this.le,
                                                bam = loosereads.bam,
                                                pad = 5000,
                                                mate_pad = 150,
                                                verbose = FALSE)
            expect_true(!is.null(this.params$ranges))
            expect_true(!is.null(this.params$qnames))
            expect_true(length(this.params$ranges) > 0)
            expect_true(length(this.params$qnames) > 0)
            expect_true(length(this.params$ranges %&% windows) > 0)
            expect_true(length(intersect(this.params$qnames, qnames)) > 0)
        }
    )
})
    

