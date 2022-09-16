{
    library(optparse)
    library(devtools)
    library(gUtils)
    library(RSeqLib)

    ## do not fail silently!
    options(error = function() {traceback(2); quit("no", 1)})

    if (!exists('opt')) {
        option_list = list(
            make_option(c("-i", "--id"), type = "character", help = "Sample / pair ID (required)"),
            make_option("--tbam", type = "character", help = "Path to tumor bam (required)"),
            make_option("--nbam", type = "character", help = "Path to normal bam (required)"),
            make_option("--ref", type = "character", help = "path to fasta for realignment (required)"),
            make_option("--concat", type = "character", help = "path to fasta for realignment to satellite and viral regions (required)", default = "/gpfs/commons/home/zchoo/git/loosends/inst/extdata/hg19_loosends/concatenated_references_deduped.fasta"),
            make_option("--loose", type = "character",
                        help = "path to loose ends file - either gGraph/JaBbA output or text file (required)"),
            make_option("--overwrite", type = "logical",
                        help = "overwrite existing analyses in output directory?",
                        default=TRUE),
            make_option("--pad", type = "numeric",
                        help = "loose end window padding",
                        default = 5e3),
            make_option("--anchorlift_window", type = "numeric",
                        help = "loose read anchorlift window",
                        default = 1e4),
            make_option("--bowtie", type = "logical", help = "use bowtie?", default=FALSE),
            make_option("--bowtie_dir", type = "character", help = "path to bowtie index directory",
                        default="~/git/loosends/inst/extdata/hg19_loosends"),
            make_option("--bowtie_ref", type = "character", help = "basename of bowtie index files",
                        default="human_g1k_v37_decoy"),
            make_option("--libdir", type = "character", help = "libdir", default=""),
            make_option("--outdir", type = "character", help = "outdir", default="./")
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
    library(gGnome)

    ## output filenames
    loose.ends.fn = paste0(opt$outdir, "/", "loose.ends.rds")
    all.reads.fn = paste0(opt$outdir, "/", "all.reads.rds")
    anchor.fn = paste0(opt$outdir, "/", "anchor.reads.rds")
    discord.fn = paste0(opt$outdir, "/", "discord.reads.rds")
    ctypes.anchors.fn = paste0(opt$outdir, "/", "ctypes.anchors.reads.rds")
    telo.anchors.fn = paste0(opt$outdir, "/", "telomere.anchors.rds")

    if (opt$overwrite | !file.exists(loose.ends.fn)) {
        ## read loose end file
        if (grepl(".rds", opt$loose, ignore.case = TRUE)) {

            loose.dt = readRDS(opt$loose)

            ## check for gGraph input
            if (inherits(loose.dt, 'list')) {
                if (all(c('segstats', 'adj', 'purity', 'ploidy') %in% names(loose.dt))) {
                    message("Detected JaBbA output!")
                    jab = gG(jabba = loose.dt)
                    loose.dt = as.data.table((jab$loose %Q% (terminal == FALSE)))
                } else {
                    stop("Invalid file supplied")
                }
            } else if (inherits(loose.dt, 'gGraph')) {
                message("Detected gGraph input!")
                loose.dt = as.data.table((loose.dt$loose %Q% (terminal == FALSE)))
            } else if (inherits(loose.dt, 'data.table')) {
                message("Detected data.table!")
            } else if (inherits(loose.dt, 'data.frame')) {
                message("Detected data.frame")
                loose.dt = as.data.table(loose.dt)
            } else if (inherits(loose.dt, 'GRanges')) {
                message("Detected GRanges!")
                loose.dt = as.data.table(loose.dt)
            } else if (inherits(loose.dt, 'GRangesList')) {
                message("Detected GRangesList!")
                loose.dt = as.data.table(unlist(loose.dt))
            } else {
                stop("Invalid file supplied for loose ends")
            }
            
        } else {
            loose.dt = fread(opt$loose)[sample == opt$id,]
        }

        ## check for nonempty table
        if (loose.dt[, .N]) {
            if (!("sample" %in% colnames(loose.dt))) {
                loose.dt[, sample := opt$id]
            }
            loose.dt = loose.dt[(seqnames %in% c(as.character(1:22), "X", "Y")) |
                                (seqnames %in% paste0('chr', c(as.character(1:22), "X", "Y"))), ]
        }

        message("Saving loose ends for this sample...")
        saveRDS(loose.dt, loose.ends.fn)
    } else {
        loose.dt = readRDS(loose.ends.fn)
    }


    ## if there are no loose ends, save empty data tables and exit
    if (!nrow(loose.dt)) {
        message("Sample has no loose ends!")
        saveRDS(data.table(), all.reads.fn)
        saveRDS(data.table(), anchor.fn)
    } else {

        if (opt$overwrite | !file.exists(all.reads.fn) | file.info(all.reads.fn)$size == 0) {
            message("Grabbing loose reads")

            if (opt$bowtie) {
                reads.dt = loosereads_wrapper(dt2gr(loose.dt),
                                              tbam = opt$tbam,
                                              nbam = opt$nbam,
                                              ref = opt$ref,
                                              id = opt$id,
                                              outdir = opt$outdir,
                                              bowtie = TRUE,
                                              bowtie.ref = opt$bowtie_ref,
                                              bowtie.dir = opt$bowtie_dir,
                                              pad = opt$pad,
                                              verbose = TRUE)
            } else {
                reads.dt = loosereads_wrapper(dt2gr(loose.dt),
                                              tbam = opt$tbam,
                                              nbam = opt$nbam,
                                              ref = opt$ref,
                                              id = opt$id,
                                              outdir = opt$outdir,
                                              pad = opt$pad,
                                              verbose = TRUE)
            }
            saveRDS(reads.dt, all.reads.fn)
        } else {
            message("Reading loose reads from file")
            reads.dt = readRDS(all.reads.fn)
            anchor.tmp.dt = reads.dt[(high.mate) & (loose.pair),]
            discord.dt = reads.dt[(!concord),]
        }

        if (opt$overwrite | !file.exists(anchor.fn) | file.info(anchor.fn)$size == 0) {
            if (!nrow(reads.dt)) {
                message("No loose reads!")
                anchor.tmp.dt = reads.dt[(high.mate) & (loose.pair),]
                discord.dt = reads.dt[(!concord),]
                saveRDS(data.table(), anchor.fn)
            } else {

                message("Annotating loose pair seeds")
                anchor.tmp.dt = reads.dt[(high.mate) & (loose.pair),]
                if (!nrow(anchor.tmp.dt)) {
                    anchor.dt = data.table()
                } else {
                    anchor.gr = dt2gr(anchor.tmp.dt[, .(seqnames, start, end, strand)],
                                      seqlengths = hg_seqlengths())
                    loose.gr = dt2gr(loose.dt[, .(seqnames, start, end, strand)],
                                     seqlengths = hg_seqlengths())

                    ## actually use gUtils anchorlift
                    anchorlift.gr = gUtils::anchorlift(query = anchor.gr,
                                                       subject = loose.gr,
                                                       window = opt$anchorlift_window)

                    anchor.dt = as.data.table(anchorlift.gr)

                    
                    ## add annotations
                    anchor.dt[, sample := anchor.tmp.dt$sample[query.id]]
                    anchor.dt[, track := anchor.tmp.dt$track[query.id]]
                    anchor.dt[, tumor := ifelse(track %like% "sample", "tumor", "normal")]
                    anchor.dt[, query.strand := anchor.tmp.dt$strand[query.id]]
                    anchor.dt[, subject.strand := loose.dt$strand[subject.id]]
                    anchor.dt[, loose.end := loose.dt$loose.end[subject.id]]
                    anchor.dt[, qname := anchor.tmp.dt$qname[query.id]]

                    ## flip the strands because the loose ends strands are in "junction orientation"
                    anchor.dt[, forward := query.strand != subject.strand]

                    ## add whether the loose end is annotated as 'germline' or 'somatic'
                    anchor.dt[, somatic := NA]
                    if (!is.null(loose.dt$keep)) {
                        anchor.dt[, somatic := ifelse(loose.dt$keep[subject.id], "somatic", "germline")]
                    }
                    
                }

                saveRDS(anchor.dt, anchor.fn)
            }
        } else {
            message("anchorlift file already exists!")
            anchor.dt = readRDS(anchor.fn)
            anchor.tmp.dt = reads.dt[(high.mate) & (loose.pair),]
            anchor.dt[, qname := anchor.tmp.dt$qname[as.numeric(as.character(query.id))]]
        }

        if (opt$overwrite | !file.exists(discord.fn) | file.info(discord.fn)$size == 0) {
            if (!nrow(reads.dt)) {
                message("No loose reads!")
                anchor.tmp.dt = reads.dt[(high.mate) & (loose.pair),]
                discord.dt = reads.dt[(!concord),]
                saveRDS(data.table(), anchor.fn)
            } else {

                message("anchorlifting discordant pairs")
                discord.dt = reads.dt[(!concord),]
                if (!nrow(discord.dt)) {
                    discord.anchorlift.dt = data.table()
                } else {
                    anchor.gr = dt2gr(discord.dt[, .(seqnames, start, end, strand)],
                                      seqlengths = hg_seqlengths())
                    loose.gr = dt2gr(loose.dt[, .(seqnames, start, end, strand)],
                                     seqlengths = hg_seqlengths())

                    ## actually use gUtils anchorlift
                    discord.anchorlift.gr = gUtils::anchorlift(query = anchor.gr,
                                                               subject = loose.gr,
                                                               window = opt$anchorlift_window)

                    discord.anchorlift.dt = as.data.table(discord.anchorlift.gr)

                    
                    ## add annotations
                    discord.anchorlift.dt[, sample := discord.dt$sample[query.id]]
                    discord.anchorlift.dt[, track := discord.dt$track[query.id]]
                    discord.anchorlift.dt[, tumor := ifelse(track %like% "sample", "tumor", "normal")]
                    discord.anchorlift.dt[, query.strand := discord.dt$strand[query.id]]
                    discord.anchorlift.dt[, subject.strand := loose.dt$strand[subject.id]]
                    discord.anchorlift.dt[, loose.end := loose.dt$loose.end[subject.id]]
                    discord.anchorlift.dt[, qname := discord.dt$qname[query.id]]

                    ## flip the strands because the loose ends strands are in "junction orientation"
                    discord.anchorlift.dt[, forward := query.strand != subject.strand]

                    ## add whether the loose end is annotated as 'germline' or 'somatic'
                    discord.anchorlift.dt[, somatic := NA]
                    if (!is.null(loose.dt$keep)) {
                        discord.anchorlift.dt[, somatic := ifelse(loose.dt$keep[subject.id], "somatic", "germline")]
                    }
                    
                }

                saveRDS(discord.anchorlift.dt, discord.fn)
            }
        } else {
            message("anchorlift file already exists!")
            discord.anchorlift.dt = readRDS(discord.fn)
            discord.dt = reads.dt[(!concord),]
            discord.anchorlift.dt[, qname := discord.dt$qname[as.numeric(as.character(query.id))]]
        }

        if (opt$overwrite | (!file.exists(ctypes.anchors.fn))) {
        ## if (TRUE) {

            loose.anchors.dt = reads.dt[(loose.pair) & (high.mate),]
            loose.mates.dt = reads.dt[(loose.pair) & (!high.mate),]

            if (!loose.mates.dt[, .N]) {
                message("No loose reads!")
                viral.anchor.dt = data.table()
                satellite.anchor.dt = data.table()
                telomeric.anchor.dt = data.table()
            } else {
                ## get all c_types of the low-mappability mate
                message("Loading reference")
                concat = RSeqLib::BWA(fasta = opt$concat)
                mate.alns.dt = as.data.table(concat[loose.mates.dt[, seq]])

                ## add back qnames
                mate.alns.dt[, qname_old := copy(qname)]
                mate.alns.dt[, qname := loose.mates.dt$qname[as.numeric(as.character(qname_old))]]

                ## replace sequence
                mate.alns.dt[, seq := loose.mates.dt$seq[as.numeric(as.character(qname_old))]]

                ## replace original strand
                mate.alns.dt[, strand := loose.mates.dt$strand[as.numeric(as.character(qname_old))]]
                
                ## check for ctypes
                mate.alns.dt = loosends::add_contig_ctypes(mate.alns.dt)

                ## add unaligned reads
                unaln.dt = loose.mates.dt[!(qname %in% mate.alns.dt$qname), ]
                mate.alns.dt = rbind(mate.alns.dt, unaln.dt, fill = TRUE)

                ## add query.seq though it really doesn't matter if the strands are flipped here
                mate.alns.dt[, query.seq := ifelse(strand == "-",
                                                   as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq))),
                                                   seq)]

                ## check for telomeres
                mate.alns.dt = check_contigs_for_telomeres(mate.alns.dt, verbose = verbose)

                ## check specifically for TTAGGG
                pattern = "TTAGGGTTAGGGTTAGGG"
                mate.alns.dt[, spec.telomeric := grepl(pattern = pattern, x = seq) |
                                   grepl(pattern = pattern, x = query.seq)]
                
                ctypes.dt = mate.alns.dt[, .(viral = any(c_type %like% "viral", na.rm = TRUE),
                                             satellite = any(c_type %like% "rep", na.rm = TRUE),
                                             g_telomeric = any(query_g_telomere, na.rm = TRUE),
                                             c_telomeric = any(query_c_telomere, na.rm = TRUE),
                                             spec.telomeric,
                                             telomeric = any(query_c_telomere | query_g_telomere, na.rm = TRUE)),
                                         by = qname]

                ## make sure that qnames are unique
                ctypes.dt = unique(ctypes.dt, by = "qname")
                anchor.dt = unique(anchor.dt, by = "qname")

                ctypes.anchors.dt = merge.data.table(anchor.dt[, .(seqnames, start, end, strand, qname,
                                                                   query.id, subject.id,
                                                                   track, tumor, query.strand, subject.strand,
                                                                   forward, somatic)],
                                                     ctypes.dt,
                                                     by = "qname")

                saveRDS(ctypes.anchors.dt, ctypes.anchors.fn)
            }
        } else {
            message("cyptes file already exists!")
        }

        ## if (opt$overwrite | (!file.exists(telo.anchors.fn))) {
        if (TRUE) {
            message("Telomere anchorlift!")

            ## checking for loose reads
            loose.reads.dt = reads.dt[(loose.pair),]

            if (loose.reads.dt[, .N]) {
                wide.reads.dt = merge.data.table(loose.reads.dt[(high.mate),
                                                          .(seqnames, start, end, strand,
                                                            mapq, track, seq, reading.frame, sample,
                                                            qname)],
                                                 loose.reads.dt[(!high.mate),
                                                          .(seqnames, start, end, strand,
                                                            mapq, track, seq, reading.frame, sample,
                                                            qname)],
                                                 by = "qname",
                                                 suffixes = c(".anchor", ".mate"))

                ## use reading frame - this was the input to the sequencer
                ## put the sequence back in the reading frame of its high mapq mate
                wide.reads.dt[is.na(reading.frame.mate), reading.frame.mate := ""]
                wide.reads.dt[, og.seq := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(reading.frame.mate)))]

                ## if the anchor strand is minus, reverse complement it
                wide.reads.dt[, query.seq := copy(og.seq)]


                ## check for G telomeres, specifically
                wide.reads.dt[, g_telomeres := loosends::find_telomeres(seq = query.seq, gorc = "g")]
                wide.reads.dt[, c_telomeres := loosends::find_telomeres(seq = query.seq, gorc = "c")]

                ## put reads in coordinates of loose ends
                alift.gr = anchorlift(subject = dt2gr(loose.dt[, .(seqnames, start, end, strand)]),
                                      query = dt2gr(wide.reads.dt[, .(seqnames = seqnames.anchor,
                                                                      start = start.anchor,
                                                                      end = end.anchor,
                                                                      strand = strand.anchor)]),
                                      window = 1e6)


                ## transfer annotations
                alift.dt = as.data.table(alift.gr)

                ## keep original loose end location
                alift.dt[, loose.end := loose.dt$loose.end[subject.id]]

                ## transfer g and c telomere annotations
                alift.dt[, g_telomeres := wide.reads.dt$g_telomeres[query.id]]
                alift.dt[, c_telomeres := wide.reads.dt$c_telomeres[query.id]]

                ## add the read
                alift.dt[, query.seq := wide.reads.dt$query.seq[query.id]]
                alift.dt[, og.seq := wide.reads.dt$og.seq[query.id]]

                ## transfer strand annotations
                alift.dt[, loose.strand := loose.dt$strand[subject.id]]
                alift.dt[, read.strand := wide.reads.dt$strand.anchor[query.id]]
                alift.dt[, mate.strand := wide.reads.dt$strand.mate[query.id]]

                ## transfer tumor/normal annotations
                alift.dt[, tumor := ifelse(wide.reads.dt$track.anchor[query.id] %like% "sample",
                                            "tumor",
                                            "normal")]

                ## transfer forward/referse annotations
                alift.dt[, forward := ifelse(loose.strand != read.strand,
                                             "forward",
                                             "reverse")]

                ## transfer germline/somatic annotations
                alift.dt[, somatic := ifelse(loose.dt$keep[subject.id],
                                             "somatic",
                                             "germline")]

                

                
            } else {
                message("No loose reads!")
                alift.dt = data.table()
            }

            saveRDS(alift.dt, telo.anchors.fn)
        } else {
            message("Telomere anchorlift complete!")
            alift.dt = readRDS(telo.anchors.fn)
        }
        
        message("Done!")
    }
    quit("no", status = 0)
}
