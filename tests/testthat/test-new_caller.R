library(loosends)

tra.bp.gr.fn = system.file("extdata", "tests", "new_caller_1", "tra.bp.gr.rds", package = "loosends")
tra.reads.dt.fn = system.file("extdata", "tests", "new_caller_1", "tra.reads.rds", package = "loosends")

tra.bp.gr = readRDS(tra.bp.gr.fn)
tra.reads.dt = readRDS(tra.reads.dt.fn)

concat.fasta.fn = system.file("extdata", "hg19_loosends", "concatenated_references_deduped.fasta", package = "loosends")
human.fasta.fn = system.file("extdata", "hg19_loosends", "human_g1k_v37_decoy.fasta", package = "loosends")
concat = BWA(fasta = concat.fasta.fn)
human = BWA(fasta = human.fasta.fn)

## TODO: STASH THIS
dj = readRDS("~/projects/gGnome/files/zc_stash/unmappable.genome.rds")

test_that(desc = "test new caller wrapper", code = {
    suppressWarnings(
        expr = {
            res = call_loose_end2(li = tra.bp.gr,
                                  ri = tra.reads.dt,
                                  id = "387",
                                  concat.bwa = concat,
                                  human.bwa = human,
                                  unmappable.gr = dj,
                                  unmappable.pad = 5e3,
                                  window = 1000,
                                  verbose = TRUE)
        })
})

