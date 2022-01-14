library(testthat)
library(loosends)
library(gUtils)

## reference directory
ref.dir = "~/git/loosends/inst/extdata/hg19_loosends"
ref.fasta = "~/git/loosends/inst/extdata/hg19_loosends/human_g1k_v37_decoy.fasta"
reads.dt.fn = system.file("extdata", "tests", "loosereads_1", "reads.dt.rds", package = "loosends")
le.dt.fn = system.file("extdata", "tests", "loosereads_1", "le.dt.rds", package = "loosends")
contig.dt.fn = system.file("extdata", "tests", "loosereads_1", "single.contig.rds", package = "loosends")
multiple.dt.fn = system.file("extdata", "tests", "loosereads_1", "multiple.contig.rds", package = "loosends")

## params
this.pair = "G32831.HCC1954"
this.le = parse.gr("20:60158837-")

## load saved data
le.dt = readRDS(le.dt.fn)
reads.dt = readRDS(reads.dt.fn)
single.dt = readRDS(contig.dt.fn)
multiple.dt = readRDS(multiple.dt.fn)

## load ref bwa
ref.obj = grab_ref_obj(ref.dir = ref.dir)

## test_that(desc = "test contig support for single contig", code = {
##     suppressWarnings(
##         expr = {
##             rc = read_support(le.dt = le.dt,
##                               reads.dt = reads.dt,
##                               contigs.dt = single.dt,
##                               ref.bwa = ref.obj$human,
##                               verbose = FALSE)
##             ## expect non-zero unmber of supporting reads
##             expect_true(rc[, .N] > 0)
##             ## expect all reads returned to be supporting reads
##             expect_true(all(rc[, supporting]))
##             expect_true(all(rc[, sample == this.pair]))
##         }
##     )
## })

## test_that(desc = "test contig support for single contig with bowtie", code = {
##     suppressWarnings(
##         expr = {
##             rc = read_support(le.dt = le.dt,
##                               reads.dt = reads.dt,
##                               contigs.dt = single.dt,
##                               fasta = ref.fasta, ##pastref.bwa = NULL,ref.obj$human,
##                               outdir = "~/testing_tmp/bowtie2_read_support",
##                               bowtie = TRUE,
##                               verbose = FALSE)
##             ## expect non-zero unmber of supporting reads
##             expect_true(rc[, .N] > 0)
##             ## expect all reads returned to be supporting reads
##             expect_true(all(rc[, supporting]))
##             expect_true(all(rc[, sample == this.pair]))
##         }
##     )
## })


## test_that(desc = "test contig support wrapper for multiple contigs", code = {
##     suppressWarnings(
##         expr = {
##             rc = read_support_wrapper(le.dt = le.dt,
##                                       reads.dt = reads.dt,
##                                       contigs.dt = multiple.dt,
##                                       id = this.pair,
##                                       ref.bwa = ref.obj$human,
##                                       verbose = FALSE)
##             expect number of rows to be equal to number of unique contigs
##             expect_true(rc[, .N] == length(unique(multiple.dt[, qname])))
##             expect all leix's to be represented
##             expect_true(all(multiple.dt[, leix] %in% rc[, leix]))
##         }
##     )
## })

test_that(desc = "test contig support wrapper with bowtie", code = {
    suppressWarnings(
        expr = {
            rc = read_support_wrapper(le.dt = le.dt,
                                      reads.dt = reads.dt,
                                      contigs.dt = multiple.dt,
                                      id = this.pair,
                                      ref.bwa = NULL,##ref.obj$human,
                                      ref.fasta = ref.fasta,
                                      bowtie = TRUE,
                                      outdir = "~/testing_tmp/contigsupport/wrapper_bowtie2",
                                      verbose = FALSE)
            ## expect number of rows to be equal to number of unique contigs
            expect_true(rc[, .N] == length(unique(multiple.dt[, qname])))
            ## expect all leix's to be represented
            expect_true(all(multiple.dt[, leix] %in% rc[, leix]))
        }
    )
})

