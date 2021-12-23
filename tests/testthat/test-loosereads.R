library(loosends)
library(testthat)

## generally useful files
ref = "~/DB/hg19/human_g1k_v37_decoy.fasta"

## TEST CASE 1
## files for first test case
loosereads.bam = system.file("tests", "testthat", "data", "loosereads_1", "loosereads.bam", package = "loosends")
aln.bam = system.file("tests", "testthat", "data", "loosereads_1", "aln.bam", package = "loosends")
normal.loosereads.bam = system.file("tests", "testthat", "data", "loosereads_1_normal", "loosereads.bam", package = "loosends")
normal.aln.bam = system.file("tests", "testthat", "data", "loosereads_1_normal", "aln.bam", package = "loosends")
qnames.txt = system.file("tests", "testthat", "data", "loosereads_1", "qnames.txt", package = "loosends")
windows.bed = system.file("tests", "testthat", "data", "loosereads_1", "windows.bed", package = "loosends")

## loci for first test case
this.le = parse.gr("20:60158837-")

## test sample name for the first test case
this.pair = "G32831.HCC1954"
this.pair.normal = "G32831.HCC1954N"

## expected results for the first test case
qnames = readLines(con = qnames.txt)
windows = rtracklayer::import.bed(con = windows.bed)

## test getting qnames + windows for reading bam
test_that(desc = "check that correct qnames and windows are found", code = {
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
    
## test bam subsetting
test_that(desc = "check that reads and mates can be recovered", code = {
    suppressWarnings(
        expr = {
            ## make sure bam file is created and has nonzero size
            subset.bam = grab_loosereads(bam = loosereads.bam,
                                         ranges = windows,
                                         qnames = qnames,
                                         outdir = "~/testing_tmp",
                                         overwrite = TRUE,
                                         verbose = FALSE)
            expect_true(file.exists(subset.bam))
            expect_true(file.info(subset.bam)$size > 0)

            ## make sure qnames and windows overlap
            subset.bam.grl = bamUtils::read.bam(subset.bam, all = TRUE,
                                               isPaired = TRUE,
                                               pairs.grl = TRUE,
                                               isDuplicate = NA)
            subset.bam.gr = unlist(subset.bam.grl)
            expect_true(length(intersect(qnames, subset.bam.gr$qname)) > 0)
            expect_true(length(subset.bam.gr %&% windows) > 0)

            ## make sure reads and mates are both there
            subset.bam.dt = as.data.table(subset.bam.gr)
            expect_true(all(subset.bam.dt[, .N, by = qname]$N == 2))
        }
    )
})

## test realignment
test_that(desc = "check realignment", code = {
    suppressWarnings(
        expr = {
            realn.bam = realign_loosereads(bam = loosereads.bam,
                                           ref = ref,
                                           outdir = "~/testing_tmp",
                                           verbose = FALSE)
            expect_true(file.exists(realn.bam))
            expect_true(file.info(realn.bam)$size > 0)

            ## check that this matches with original realignment
            aln.bam.grl = bamUtils::read.bam(aln.bam, all = TRUE,
                                             isPaired = NA,
                                             pairs.grl = TRUE,
                                             isDuplicate = NA)
            aln.bam.gr = unlist(aln.bam.grl)
            realn.bam.grl = bamUtils::read.bam(realn.bam, all = TRUE,
                                                isPaired = NA,
                                                pairs.grl = TRUE,
                                                isDuplicate = NA)
            realn.bam.gr = unlist(realn.bam.grl)
            expect_true(length(intersect(realn.bam.gr$qname, aln.bam.gr$qname)) > 0)
            expect_true(length(realn.bam.gr %&% aln.bam.gr) > 0)
        }
    )
})

## test loose read merging and annotation
test_that(desc = "check loose read merging and annotation", code = {
    suppressWarnings(
        expr = {
            loosereads.res = loose.reads2(tbam = loosereads.bam,
                                          taln = aln.bam,
                                          filter = FALSE,
                                          verbose = FALSE)

            ## check that expected qnames are there
            expect_true(length(intersect(qnames, loosereads.res$qname)) > 0)
            ## check that read and mate are both present
            expect_true(all(loosereads.res[, .N, by = qname]$N == 2))
            
            filtered.res = loose.reads2(tbam = loosereads.bam,
                                        taln = aln.bam,
                                        filter = TRUE,
                                        verbose = FALSE)

            ## check that expected qnames are there
            expect_true(length(intersect(qnames, filtered.res$qname)) > 0)
            ## check that read and mate are both present
            expect_true(all(filtered.res[, .N, by = qname]$N == 2))
            ## check that all reads returned are loose
            expect_true(all(filtered.res[, loose.pair]))
            expect_true(all(filtered.res[(high.mate), mapq] >= 50))
            expect_true(all(filtered.res[(!high.mate), is.na(mapq) | mapq == 0]))
        }
    )
})

## test loose read merging and annotation with matched normal
test_that(desc = "check loose read annotation with matched normal", code = {
    suppressWarnings(
        expr = {
            full.loosereads.res = loose.reads2(tbam = loosereads.bam,
                                          taln = aln.bam,
                                          nbam = normal.loosereads.bam,
                                          naln = normal.aln.bam,
                                          id = this.pair,
                                          filter = FALSE,
                                          verbose = FALSE)

            ## check that mate and seed are present
            expect_true(all(full.loosereads.res[, .N, by = .(sample, qname)]$N == 2))
            ## check that sample and matched normal both present
            expect_true(this.pair %in% full.loosereads.res$sample &
                        this.pair.normal %in% full.loosereads.res$sample)
        }
    )
})
