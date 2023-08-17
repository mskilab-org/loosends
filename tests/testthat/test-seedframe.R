library(loosends)

inv.bp.gr.fn = system.file("extdata", "tests", "new_caller_1", "inv.bp.gr.rds", package = "loosends")
inv.reads.dt.fn = system.file("extdata", "tests", "new_caller_1", "inv.reads.rds", package = "loosends")

inv.bp.gr = readRDS(inv.bp.gr.fn)
inv.reads.dt = readRDS(inv.reads.dt.fn)

test_that(desc = "check seed frame preparation for an FBI", code = {
    suppressWarnings(
        expr = {
            le.dt = prep_loose_ends(li = inv.bp.gr, id = "2527")
            prepped.reads.dt = prep_loose_reads(li = le.dt, loose.reads.dt = inv.reads.dt)
            seed.rds = grab_seed_frame(prepped.reads.dt,
                                       gr.flipstrand(inv.bp.gr) + 1e3,
                                       seq.field = "reading.frame")
            ## check that there is one seed read per qname pair
            expect_true(all(seed.rds[, sum(seed), by = qname]$V1 <= 1))
            ## check that seed reads have the correct strand
            expect_true(all(seed.rds[(seed), as.character(strand)] == strand(gr.flipstrand(inv.bp.gr))))
            ## check that seed reads are not reverse complemented, but non-seed reads are
            expect_true(all(seed.rds[(seed), reading.frame == seed.frame]))
            expect_true(all(seed.rds[(!seed), reading.frame != seed.frame]))
        })
})

dup.bp.gr.fn = system.file("extdata", "tests", "new_caller_1", "dup.bp.gr.rds", package = "loosends")
dup.reads.dt.fn = system.file("extdata", "tests", "new_caller_1", "dup.reads.rds", package = "loosends")

dup.bp.gr = readRDS(dup.bp.gr.fn)
dup.reads.dt = readRDS(dup.reads.dt.fn)

test_that(desc = "check seed frame preparation for a DUP", code = {
    suppressWarnings(
        expr = {
            le.dt = prep_loose_ends(li = dup.bp.gr, id = "1109")
            prepped.reads.dt = prep_loose_reads(li = le.dt, loose.reads.dt = dup.reads.dt)
            seed.rds = grab_seed_frame(prepped.reads.dt,
                                       gr.flipstrand(dup.bp.gr) + 1e3,
                                       seq.field = "reading.frame")
            ## check that there is one seed read per qname pair
            expect_true(all(seed.rds[, sum(seed), by = qname]$V1 <= 1))
            ## check that seed reads have the correct strand
            expect_true(all(seed.rds[(seed), as.character(strand)] == strand(gr.flipstrand(dup.bp.gr))))
            ## check that seed reads are not reverse complemented, but non-seed reads are
            expect_true(all(seed.rds[(seed), reading.frame == seed.frame]))
            expect_true(all(seed.rds[(!seed), reading.frame != seed.frame]))
        })
})

tra.bp.gr.fn = system.file("extdata", "tests", "new_caller_1", "tra.bp.gr.rds", package = "loosends")
tra.reads.dt.fn = system.file("extdata", "tests", "new_caller_1", "tra.reads.rds", package = "loosends")

tra.bp.gr = readRDS(tra.bp.gr.fn)
tra.reads.dt = readRDS(tra.reads.dt.fn)

test_that(desc = "check seed frame preparation for a TRA", code = {
    suppressWarnings(
        expr = {
            le.dt = prep_loose_ends(li = tra.bp.gr, id = "387")
            prepped.reads.dt = prep_loose_reads(li = le.dt, loose.reads.dt = tra.reads.dt)
            seed.rds = grab_seed_frame(prepped.reads.dt,
                                       gr.flipstrand(tra.bp.gr) + 1e3,
                                       seq.field = "reading.frame")
            ## check that there is one seed read per qname pair
            expect_true(all(seed.rds[, sum(seed), by = qname]$V1 <= 1))
            ## check that seed reads have the correct strand
            expect_true(all(seed.rds[(seed), as.character(strand)] == strand(gr.flipstrand(tra.bp.gr))))
            ## check that seed reads are not reverse complemented, but non-seed reads are
            expect_true(all(seed.rds[(seed), reading.frame == seed.frame]))
            expect_true(all(seed.rds[(!seed), reading.frame != seed.frame]))
        })
})

          
