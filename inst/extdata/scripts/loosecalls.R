{
    library(optparse)
    library(devtools)
    
    ## DO NOT FAIL SILENTLY
    options(error = function() {traceback(2); quit("no", 1)})

    if (!exists('opt')) {
        option_list = list(
            make_option(c("-i", "--id"), type = "character", help = "Sample / pair ID (required)"),
            make_option("--libdir", type = "character", help = "libdir", default=""),
            make_option("--outdir", type = "character", help = "outdir", default="./"),
            make_option("--ref_dir", type = "character", help = "ref dir", default="./"),
            make_option("--pad", type = "numeric", help = "assembly window (bp)", default=5e3),
            make_option("--max.qnames", type = "numeric", help = "max number of qnames", default=1e4),
            make_option("--mix", type = "logical", help = "mix t/n reads before assembly", default=FALSE),
            make_option("--overwrite", type = "logical", help = "overwrite existing analysis", default=FALSE),
            make_option("--minimap", type = "logical", help = "use minimap?", default=FALSE),
            make_option("--maxwin", type = "numeric", help = "max window around loose end to grab reads", default=5e3),
            make_option("--reads", type = "character", help = "data table of reads", default = "/dev/null"),
            make_option("--recall", type = "logical", help = "re-call?", default = TRUE),
            make_option("--cores", type = "numeric", help = "number of cores", default = 8),
            make_option("--loose", type = "character", help = "data table of loose ends", default="/dev/null")
        )
        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)

        print(opt)
        
        print(.libPaths())
        options(error=function()traceback(2))

        ## keep record of run
        writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }

    message("Loading development version of loosends")
    devtools::load_all("~/git/loosends")
    library(gUtils)

    ## output filenames
    calls.dt.fn = file.path(opt$outdir, "calls.rds")
    recall.dt.fn = file.path(opt$outdir, "recalls.rds")
    contigs.dt.fn = file.path(opt$outdir, "contigs.rds")
    filtered.dt.fn = file.path(opt$outdir, "filtered.rds")


    ## read loose end file
    if (grepl(".rds", opt$loose, ignore.case = TRUE)) {
        loose.dt = readRDS(opt$loose)
        if (loose.dt[, .N]) {
            if (!("sample" %in% colnames(loose.dt))) {
                loose.dt[, sample := opt$id]
            }
            loose.dt = loose.dt[seqnames %in% c(as.character(1:22), "X", "Y"),]
        }
    } else {
        loose.dt = fread(opt$loose)[sample == opt$id,]
    }

    if (!nrow(loose.dt)) {
        message("No loose ends for this sample!")

        call.res = data.table()
        filtered.res = data.table()
        contig.res = data.table()
    } else {
        ## add required metadata
        message("Prepping loose ends")
        ## ALERT!!!!
        loose.dt = loose.dt[, sample := opt$id]
        loose.dt = prep_loose_ends(li = loose.dt, id = opt$id)

        ## grab loose reads
        message("Grabbing reads")
        if (grepl(".rds", opt$reads, ignore.case = TRUE)) {
            loose.reads.dt = readRDS(opt$reads)
        } else {
            loose.reads.dt = fread(opt$reads)
        }

        ## convert reads to GRanges
        if (inherits(loose.reads.dt, 'data.table')) {
            loose.reads.gr = dt2gr(loose.reads.dt[, .(seqnames, start, end)])
            loose.ends.gr = dt2gr(loose.dt[, .(seqnames, start, end, sample, leix)])
        } else {
            stop("loose reads must be data.table")
        }

        ## load reference from disk
        message("Loading references from disk")

        ## only run if filtered contigs don't already exist
        if (opt$overwrite || (!file.exists(filtered.dt.fn))) {

            human.ref = paste0(opt$ref_dir, "/human_g1k_v37_decoy.fasta")
            concat.ref = paste0(opt$ref_dir, "/concatenated_references_deduped.fasta")

            human.bwa = RSeqLib::BWA(human.ref)
            concat.bwa = RSeqLib::BWA(concat.ref)

            ## if running minimap, pre-index the fasta to reduce some overhead
            if (opt$minimap) {
                message("Indexing minimap reference!")
                mmi = minimap_index(fasta = concat.ref, outdir = opt$outdir, verbose = TRUE)
            } else {
                message("Skipping indexing! Not using minimap.")
                mmi = concat.ref
            }

            message("Starting caller")
            res = mclapply(1:loose.dt[, .N],
                           function(ix, reads, loose) {
                               message("Starting analysis for loose end ", ix,
                                       " of ",
                                       loose.dt[, .N])
                               distances = GenomicRanges::distance(loose.ends.gr[ix],
                                                                   loose.reads.gr,
                                                                   ignore.strand = TRUE)
                               qnames = reads[which(distances < opt$maxwin),qname]
                               n.qnames = length(unique(qnames))
                               message("Number of unique qnames: ", n.qnames)
                               if (n.qnames > opt$max.qnames) {
                                   message("Downsampling to ", opt$max.qnames, " qnames")
                                   qnames = sample(qnames, size = opt$max.qnames, replace = FALSE)
                               }
                               message("Prepping reads...")
                               this.loose.reads = reads[qname %in% qnames,]
                               ri = prep_loose_reads(loose[ix,], this.loose.reads)
                               message("Starting caller!")
                               call.res = loosends:::call_loose_end2(loose[ix,],
                                                                     ri,
                                                                     id = opt$id,
                                                                     concat.bwa = concat.bwa,
                                                                     human.bwa = human.bwa,
                                                                     window = opt$pad,
                                                                     verbose = (opt$cores == 1))
                               message("Annotation: ", call.res$label.dt$annotation)
                               message("Somatic: ", call.res$label.dt$somatic)
                               return(call.res)
                           },
                           loose.reads.dt,
                           loose.dt,
                           mc.cores = opt$cores)

            message("Concatenating calls")
            call.res = lapply(res, function(x) {return(x$label.dt)}) %>%
                rbindlist(fill = TRUE)
            contig.res = lapply(res, function(x) {return(x$all.contigs)}) %>%
                rbindlist(fill = TRUE)
            filtered.res = lapply(res, function(x) {return(x$keep.tigs.support)}) %>%
                rbindlist(fill = TRUE)

            call.res[, sample := opt$id]

            if (filtered.res[, .N]) {
                filtered.res[, sample := opt$id]
            }
            message("Saving output")
            saveRDS(call.res, calls.dt.fn)
            saveRDS(contig.res, contigs.dt.fn)
            saveRDS(filtered.res, filtered.dt.fn)
        } else {
            message("Filtered contigs already exist!!")
            filtered.res = readRDS(filtered.dt.fn)
        }
    }

    

    ## check if filtered contigs should be reprocessed
    if (opt$overwrite || opt$recall || !file.exists(recall.dt.fn)) {
        if (loose.dt[, .N]) {
            reres = lapply(1:loose.dt[, .N],
                           function(ix, filtered, loose) {
                               message("Recalling loose end: ", ix)
                               ## do this based on seed overlaps
                               if (!is.null(filtered$loose.end)) {
                                   calns = filtered[(loose.end == loose$loose.end[ix]),]
                               } else {
                                   sel = which(parse.gr(filtered[, seed]) %^^%
                                               (gr.flipstrand(parse.gr(loose[ix, loose.end]) + 5e3)))
                                   filtered[sel, loose.end := loose$loose.end[ix]]
                                   calns = filtered[sel,]
                               }
                               new.call = loosends::check_contig_concordance(calns = calns)
                               new.call[, loose.end := loose$loose.end[ix]]
                               return(new.call)
                           },
                           filtered.res,
                           loose.dt)
            reres = rbindlist(reres, fill = TRUE)
            reres[, sample := opt$id]
        } else {
            reres = data.table()
        }

        message("Saving output!")
        saveRDS(reres, recall.dt.fn)
        saveRDS(filtered.res, filtered.dt.fn)
    } else {
        message("Already re-called!")
    }

    message("Done")

    quit("no")
}
