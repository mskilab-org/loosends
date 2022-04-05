#' @name build_contigs_wrapper
#' @title build_contigs_wrapper
#'
#' @param gr (GRanges) which breakpoint do we care about
#' @param reads.dt (data.table)
#' @param ref (BWA object)
#' @param window (numeric) how many bases around window to consider for assembly (default 5e3)
#' @param assembly.region (numeric) how many bases at a time to consider for assembly (default 1e3)
#' @param stride (numeric) stride for tiling window (default 500)
#' @param pseudo.contigs (logical) include pseudo-contigs from discordant read clustering? default TRUE
#' @param low.mappability.gr (GRanges) path to GRanges containing low mappability bases
#' @param unassembled.gr (GRanges) path to GRanges containing unassembled regions
#' @param verbose (logical) default FALSE
#'
#' @return data.table with contig alignments
build_contigs_wrapper = function(gr, reads.dt, ref,
                                 window = 5e3,
                                 assembly.region = 1e3,
                                 stride = 500,
                                 pseudo.contigs = TRUE,
                                 low.mappability.gr = GRanges(),
                                 unassembled.gr = GRanges(),
                                 verbose = FALSE) {

    tiles = unlist(GenomicRanges::slidingWindows(x = gr + window,
                                                 width = assembly.region,
                                                 step = stride))

    all.contigs = lapply(seq_along(tiles),
                  function(ix) {
                      if (verbose) {message("Starting analysis for window: ", ix, " of ", length(tiles))}
                      ## get reads for forward and reverse frame separately
                      forward.seed.frame.dt = grab_seed_frame(reads.dt,
                                                              seed.gr = tiles[ix],
                                                              seq.field = "reading.frame",
                                                              forward = TRUE)
                      reverse.seed.frame.dt = grab_seed_frame(reads.dt,
                                                              seed.gr = gr.flipstrand(tiles[ix]),
                                                              seq.field = "reading.frame",
                                                              forward = FALSE)
                      if (forward.seed.frame.dt[, .N]) {
                          ## if (verbose) { message("Building forward track contigs") }
                          forward.ctigs = build_contigs(forward.seed.frame.dt, verbose = verbose)
                      } else {
                          forward.ctigs = character()
                      }

                      ctigs = forward.ctigs

                      ## only use forward strand
                      ## if (reverse.seed.frame.dt[, .N]) {
                      ##     ## if (verbose) { message("Building reverse track contigs") }
                      ##     reverse.ctigs = build_contigs(reverse.seed.frame.dt, verbose = verbose)
                      ## } else {
                      ##     reverse.ctigs = character()
                      ## }
                      ## ctigs = c(forward.ctigs, reverse.ctigs)
                          
                      if (length(ctigs)) {
                          ## if (verbose) { message("Aligning contigs to reference") }
                          aln.ctigs = align_contigs(ctigs,
                                                    ref,
                                                    verbose = verbose,
                                                    keep.unaligned = FALSE)
                          aln.ctigs = add_contig_ctypes(aln.ctigs, verbose = verbose)
                          aln.ctigs = check_contigs_for_telomeres(aln.ctigs, verbose = verbose)
                          ## if (verbose) { message("Filtering contigs based on structure") }
                          qc.ctigs = qc_contigs(aln.ctigs,
                                                tiles[ix],
                                                low.mappability.gr = low.mappability.gr,
                                                unassembled.gr = unassembled.gr)
                          if (qc.ctigs[, .N]) {
                              qc.ctigs[, name := paste(ix, qname, sep = ".")]
                              qc.ctigs[, seed := gr.string(tiles[ix])]
                          }
                      } else {
                          qc.ctigs = data.table()
                      }

                      if (pseudo.contigs) {
                          ## if (verbose) {
                          ##     message("Building pseudo-contigs from discordant read pairs")
                          ## }
                          forward.pseudo.ctigs = build_pseudo_contigs(forward.seed.frame.dt,
                                                                      verbose = verbose)
                          pseudo.ctigs = forward.pseudo.ctigs

                          ## only use forward strand reads
                          ## reverse.pseudo.ctigs = build_pseudo_contigs(reverse.seed.frame.dt,
                          ##                                             verbose = verbose)
                          ## pseudo.ctigs = c(forward.pseudo.ctigs, reverse.pseudo.ctigs)
                          
                          if (length(pseudo.ctigs)) {
                              ## if (verbose) { message("Aligning pseudo-contigs to reference") }
                              aln.pseudo.ctigs = align_contigs(pseudo.ctigs,
                                                               ref,
                                                               verbose = verbose,
                                                               keep.unaligned = FALSE)
                              aln.pseudo.ctigs = add_contig_ctypes(aln.pseudo.ctigs, verbose = verbose)
                              aln.pseudo.ctigs = check_contigs_for_telomeres(aln.pseudo.ctigs, verbose = verbose)
                              ## if (verbose) {message("Filtering pseudo-contigs by structure")}
                              qc.pseudo.ctigs = qc_contigs(aln.pseudo.ctigs,
                                                           tiles[ix],
                                                           low.mappability.gr = low.mappability.gr,
                                                           unassembled.gr = unassembled.gr)
                              if (qc.pseudo.ctigs[, .N]) {
                                  qc.pseudo.ctigs[, name := paste("pseudo", ix, qname, sep = ".")]
                                  qc.pseudo.ctigs[, seed := gr.string(tiles[ix])]
                              }
                          } else {
                              qc.pseudo.ctigs = data.table()
                          }
                          ## if (verbose) {
                          ##     message("Number of pseudo-contigs: ", qc.pseudo.ctigs[, .N])
                          ## }
                          qc.ctigs = rbind(qc.ctigs, qc.pseudo.ctigs, fill = TRUE)
                      }
                      return(qc.ctigs)
                  })

    all.contigs = rbindlist(all.contigs, fill = TRUE, use.names = TRUE)
    return(all.contigs)
}

#' @name check_contigs_for_telomeres
#' @title check_contigs_for_telomeres
#'
#' @description
#'
#' check contigs for telomeres
#'
#' @param calns (data.table from contig alignments)
#' @param verbose (logical) default FALSE
#'
#' @return data.table with columns:
#' - query_c_telomere
#' - query_g_telomere
#' - aln_c_telomere
#' - aln_g_telomere
#' - c_telomere (any alignment for the qname contains a C telomere)
#' - g_telomere (any alignment for the qname contains a G telomere)
check_contigs_for_telomeres = function(calns, verbose = FALSE) {
    calns = copy(calns)

    if (!calns[, .N]) {
        calns[, ":="(query_c_telomere = logical(),
                     query_g_telomere = logical(),
                     aln_c_telomere = logical(),
                     aln_g_telomere = logical(),
                     c_telomere = logical(),
                     g_telomere = logical())]
        return(calns)
    }
                     
    calns[, query_c_telomere := find_telomeres(seq = query.seq, gorc = "c")]
    calns[, query_g_telomere := find_telomeres(seq = query.seq, gorc = "g")]
    calns[, aln_c_telomere := find_telomeres(seq = seq, gorc = "c")]
    calns[, aln_g_telomere := find_telomeres(seq = seq, gorc = "g")]
    calns[, c_telomere := any(query_c_telomere, na.rm = TRUE), by = qname]
    calns[, g_telomere := any(query_g_telomere, na.rm = TRUE), by = qname]

    return(calns)
    
}

#' @name add_contig_ctypes
#' @title add_contig_ctypes
#'
#' @description
#'
#' Given data.table of contig alignments, add c_type and c_spec labels based on reference sequence
#'
#' @param calns (data.table of contig alignments)
#' @param refseq.fn (path to key of reference sequences)
#' @param remove.decoy (logical) remove alignments to decoy sequences? default FALSE
#' @param verbose (logical) default FALSE
#'
#' @return data.table with new columns c_type and c_spec
add_contig_ctypes = function(calns,
                             refseq.fn = system.file("extdata", "reference.sequences.rds", package = "loosends"),
                             remove.decoy = FALSE,
                             verbose = FALSE) {

    calns = copy(calns)
    
    if (!calns[, .N]) {
        calns[, c_type := character()]
        calns[, c_spec := character()]
        return(calns)
    }
    
    refseq.dt = readRDS(refseq.fn)

    calns[, c_type := refseq.dt[as.character(calns$seqnames), c_type]]
    calns[, c_spec := refseq.dt[as.character(calns$seqnames), c_spec]]

    if (remove.decoy) {
        if (verbose) { message("Checking for alignments to decoy") }
        decoy.qnames = calns[, .(has.decoy = any((c_type %like% "decoy" | c_spec %like% "decoy"))), by = qname]
        if (decoy.qnames[(has.decoy), .N]) {
            calns = calns[!(qname %in% decoy.qnames[(has.decoy), qname]),]
        }
    }

    return(calns)
}
                             
#' @name align_contigs
#' @title align_contigs
#'
#' @description
#'
#' Align contigs (character vector) to reference
#' Then annotate the segments of the alignment according to reference seqnames
#'
#' Relies on seqnames --> annotation mapping included as part of the package...
#'
#' @param tigs (character)
#' @param ref (RSeqLib::BWA)
#' @param refseq.fn (character) path to seqnames annotation file (included)
#' @param primary.only (logical) default TRUE
#' @param keep.unaligned (logical) keep unaligend contigs?? default FALSE
#' @param verbose (logical) default FALSE
#'
#' @return data.table with contig alignment and metadata fields seq, query.seq
align_contigs = function(tigs, ref,
                         primary.only = TRUE,
                         remove.decoy = TRUE,
                         keep.unaligned = FALSE,
                         verbose = FALSE) {

    if (verbose) { message("Finding reference alignments") }
    
    aln.tigs = ref[tigs]
    
    if (length(aln.tigs)) {
        mcols(aln.tigs)[, "query.seq"] = tigs[as.numeric(mcols(aln.tigs)[, "qname"])]
    }
    ##     mcols(aln.tigs)[, "c_type"] = refseq.dt[as.character(seqnames(aln.tigs)), c_type]
    ##     mcols(aln.tigs)[, "c_spec"] = refseq.dt[as.character(seqnames(aln.tigs)), c_spec]
    ##     

    ##     if (verbose) { message("Searching alignments for C and G telomeres") }
    ##     ## ## does the aligned chunk specifically have a C/G telomere?
    ##     mcols(aln.tigs)[, "aln_c_telomere"] = find_telomeres(seq = mcols(aln.tigs)[, "seq"], gorc = "c")
    ##     mcols(aln.tigs)[, "aln_g_telomere"] = find_telomeres(seq = mcols(aln.tigs)[, "seq"], gorc = "g")
    ##     ## ## does the query overall have a C/G telomere?
    ##     ## mcols(aln.tigs)[, "query_c_telomere"] = find_telomeres(seq = mcols(aln.tigs)[, "query.seq"], gorc = "c")
    ##     ## mcols(aln.tigs)[, "query_g_telomere"] = find_telomeres(seq = mcols(aln.tigs)[, "query.seq"], gorc = "g")

    ## convert to data table
    aln.tigs.dt = as.data.table(aln.tigs)

    ## merge with non-aligned contigs
    if (keep.unaligned) {
        if (verbose) {message("Checking for unaligned contigs")}
        unaln.tigs = setdiff(1:length(tigs), aln.tigs.dt[, as.numeric(qname)])
        if (length(unaln.tigs)) {
            unaln.tigs.dt = data.table(seqnames = NA_character_,
                                       start = 1,
                                       end = 0,
                                       cigar = NA_character_,
                                       seq = NA_character_,
                                       query.seq = tigs[unaln.tigs],
                                       qwidth = nchar(tigs[unaln.tigs]),
                                       qname = unaln.tigs)
            aln.tigs.dt = rbind(aln.tigs.dt, unaln.tigs.dt, use.names = TRUE, fill = TRUE)
        }
        if (verbose) {message("Number of unaligned contigs: ", length(unaln.tigs))}
    }

    ## if (aln.tigs.dt[, .N]) {

    ##     if (verbose) { message("Searching query sequences for C and G telomeres") }
        
    ##     aln.tigs.dt[, query_c_telomere := find_telomeres(seq = mcols(aln.tigs)[, "query.seq"], gorc = "c")]
    ##     aln.tigs.dt[, query_g_telomere := find_telomeres(seq = mcols(aln.tigs)[, "query.seq"], gorc = "g")]
    ## }

    if (aln.tigs.dt[, .N]) {
        if (verbose) { message("Annotating primary alignments") }
        aln.tigs.dt[, primary := bamUtils::bamflag(flag)[, "isNotPrimaryRead"] == 0]
        if (primary.only) {
            aln.tigs.dt = aln.tigs.dt[(primary),]
        }
        ## if (verbose) {
        ##     message("Checking for alignments to decoy sequences")
        ## }
        ## aln.tigs.dt[, decoy := as.character(seqnames) %in% refseq.dt[c_type == "decoy", seqnames]]
        ## if (remove.decoy) {
        ##     aln.tigs.dt = aln.tigs.dt[(!decoy)]
        ## }
    }

    return(aln.tigs.dt)
}

#' @name grab_seed_frame
#' @title grab_seed_frame
#'
#' @description
#'
#' Get reads and their mates overlapping the seed window
#' Designate one member of the pair the "seed" read and the other pair the "mate" read
#' Reverse complement the "mate" read so that the reference frame is the same
#'
#' @param reads.dt (data.table) must have columns seqnames, start, end, strand, reading.frame, qname
#' @param seed.gr (GRanges) (stranded) seed window for assembly
#' @param seq.field (character) where is the seq found? default reading.frame
#' @param forward (logical) forward track? if FALSE, will reverse complement the seed frame (default TRUE)
#' @param max.n (numeric) maximum number of low-quality "N" bases (default 0)
#' @param verbose (logical) default FALSE
grab_seed_frame = function(reads.dt, seed.gr,
                           seq.field = "reading.frame",
                           forward = TRUE,
                           max.n = 0,
                           verbose = FALSE) {

    if (!reads.dt[, .N]) {
        return(reads.dt)
    }

    ## grab qnames overlapping target window
    seed.reads = dt2gr(reads.dt[, .(seqnames, start, end, strand)]) %^^% seed.gr
    ## seed.reads = dt2gr(reads.dt[, .(seqnames, start, end, strand)]) %^% seed.gr
    seed.qnames = reads.dt[which(seed.reads), qname]
    valid.reads.dt = reads.dt[qname %in% seed.qnames,]
    valid.reads.dt[, seed := dt2gr(valid.reads.dt) %^^% seed.gr]
    ## valid.reads.dt[, seed := dt2gr(valid.reads.dt) %^% seed.gr]

    ## arbitrarily choose a read as anchor if the qname is duplicated
    ## this occurs for foldback-inversions
    ## valid.reads.dt[(seed), .N, by = qname][order(N)]
    valid.reads.dt[(seed), seed := 1:.N == 1, by = qname]
    ## stranded.seed indicates whether the seed read of each pair is on the same strand as the seed
    seed.strand = unique(as.character(strand(seed.gr)))[1]
    valid.reads.dt[, stranded.seed := strand[.SD$seed] == seed.strand, by = qname]

    ## if assembling forward track
    ## for non-seed reads, reverse complement
    ## if assembling reverse track
    ## then RC seed reads, but leave the mates as is
    valid.reads.dt[, rc.reading.frame := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(get(seq.field))))]
    if (forward) {
        valid.reads.dt[, seed.frame := ifelse(seed, get(seq.field), rc.reading.frame)]
    } else {
        valid.reads.dt[, seed.frame := ifelse(seed, rc.reading.frame, get(seq.field))]
    }

    ## get rid of anything with too many N's or zero length strings
    valid.reads.dt = valid.reads.dt[stringr::str_count(get(seq.field), "N") <= max.n & nchar(get(seq.field)) > 0,]

    return(valid.reads.dt)
}

#' @name build_pseudo_contigs
#' @title build_pseudo_contigs
#'
#' @description
#'
#' Check whether there are three (default) read pairs where the mate aligns within 1 kbp of each other
#'
#' If so returns the coordinates of the MATE (mate only)
#'
#' @param reads.dt (must have column seed frame or $col containing sequence and column seed or $scol containing whether or not the read is a seed read)
#' @param col (character) field with character of seed reading frame, default "seed.frame"
#' @param qcol (character) field indicating base quality, default "qual"
#' @param scol (character) field indicating whether read is seed, default "seed"
#' @param ccol (character) field indicating whether read is concordant, default "concord"
#' @param min.mapq (numeric) minimum mapq (default 1)
#' @param pad (numeric) pad for reducing alignments, default 1e3
#' @param min.count.thresh (numeric) minimum number of reads aligning to the pseudo-contig (default 3)
#' @param verbose (logical) default FALSE
#'
#' @return character vector of candidate contig sequences
build_pseudo_contigs = function(reads.dt,
                                col = "seed.frame",
                                qcol = "qual",
                                scol = "seed",
                                ccol = "concord",
                                min.mapq = 30,
                                pad = 1e3,
                                min.count.thresh = 5,
                                verbose = FALSE) {

    if (!reads.dt[, .N]) {
        return(character())
    }

    if (!reads.dt[!(get(scol)), .N]) {
        return(character())
    }

    if (!reads.dt[!(get(scol)),][!(get(ccol)),.N]) {
        return(character())
    }

    if (!reads.dt[!(get(scol)),][!(get(ccol)),][mapq > min.mapq, .N]) {
        return(character())
    }

    ## non-concordant reads
    reads.dt = reads.dt[!(get(ccol)),][mapq > min.mapq,]
    
    ## grab mate ranges
    ## strands are flipped because these are non-seed reads so we want to flip them to
    ## the same reading frame as the seed
    mate.gr = gr.flipstrand(dt2gr(reads.dt[!(get(scol)),][!(get(ccol)),]))
    seed.gr = gr.flipstrand(dt2gr(reads.dt[(get(scol)),][!(get(ccol)),]))
    mate.ranges = GenomicRanges::reduce(gUtils::gr.sum(mate.gr) %Q% (score > 0))

    ## check overlap with non-seed mate, and no overlap with seed mate
    if (length(mate.ranges)) {
        mate.ranges = mate.ranges[mate.ranges %NN% (gr.flipstrand(mate.gr) + pad) > min.count.thresh & mate.ranges %NN% (seed.gr + pad) == 0]
    } else {
        return(character())
    }

    if (length(mate.ranges)) {
        ## get only the ranges that are in standard sequences
        mate.ranges = dt2gr(as.data.table(mate.ranges)[grepl('^(chr)*[0-9XY]+$', as.character(seqnames)),])
        seqs = Biostrings::getSeq(Hsapiens, trim(gr.fix(gr.chr(mate.ranges), Hsapiens)))

        ## make sure that there are a decent number of read alignments to the contig
        ## and flip things where more reads align to the negative than to the positive strand
        ## of the contig
        tigs.bwa = BWA(seq = seqs)
        ralns = tigs.bwa[reads.dt[, get(col)]]

        ## this is the same logic as build_contigs
        realns.strand.dt = as.data.table(ralns)[, .(n.positive.alns = sum(strand == "+"),
                                                    n.negative.alns = sum(strand == "-")),
                                                by = seqnames]
        rc.seqnames = realns.strand.dt[n.negative.alns > n.positive.alns,
                                       as.numeric(as.character(seqnames))]
        good.seqnames = realns.strand.dt[pmax(n.positive.alns, n.negative.alns) >
                                         20 * pmin(n.positive.alns, n.negative.alns),
                                         as.numeric(as.character(seqnames))]
        seqs[rc.seqnames] = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs[rc.seqnames])))
        seqs = seqs[good.seqnames]
        return(seqs)
    } else {
        return(character())
    }
}

#' @name build_contigs
#' @title build_contigs
#'
#' @description
#'
#' Assemble reads surrounding seed window supplied in seed.frame column of data table
#' Performs alignment for a max number of iterations or until all reads are aligned
#'
#' @param reads.dt (must have column seed.frame)
#' @param col (column name containing sequences for alignment) e.g. seed.frame
#' @param qcol (column name for base quality, default 'qual')
#' @param max.iter (numeric) max number of iterations for assembly
#' @param max.unaln (numeric) max number of unaligned reads after assembly
#' @param verbose (logical) default FALSE
build_contigs = function(reads.dt, col = "seed.frame", qcol = "qual", max.iter = 3, max.unaln = 3, verbose = FALSE) {

    ## get rid of NA reads
    reads.dt = reads.dt[!is.na(get(col)),]

    if (!reads.dt[, .N]) {
        return(character())
    }


    ## TRY THREE CONDITIONS
    ## assembly everything
    ## if (verbose) {
    ##     message("Trying assembly with all reads")
    ## }
    tigs = RSeqLib::Fermi(reads = reads.dt[, get(col)], qual = reads.dt[, get(qcol)], assemble = TRUE)
    all.tigs = RSeqLib::contigs(tigs)

    ## assemble discordant pairs only
    if (reads.dt[(!concord), .N] > 4) {
        ## if (verbose) {message("Trying assembly with only discordant reads")}
        discordant.tigs = RSeqLib::Fermi(reads = reads.dt[(!concord), get(col)],
                                         qual = reads.dt[(!concord), get(qcol)],
                                         assemble = TRUE)
        all.tigs = c(all.tigs, RSeqLib::contigs(discordant.tigs))
    }

    ## assemble loose pairs only
    if (reads.dt[(loose.pair), .N] > 4) {
        ## browser()
        ## if (verbose) {message("Trying assembly with only loose reads")}
        loose.tigs = RSeqLib::Fermi(reads = reads.dt[(loose.pair), get(col)],
                                    qual = reads.dt[(loose.pair), get(qcol)],
                                    assemble = TRUE)
        all.tigs = c(all.tigs, RSeqLib::contigs(loose.tigs))
    }

    ## then align reads to all of the contigs found so far
    ## if there are unaligned reads, retry assembly just from the unaligned reads
    if (length(all.tigs)) {
        ## if (verbose) { message("Checking for reads that don't align to any contigs") }
        tigs.bwa = BWA(seq = all.tigs)
        ralns = tigs.bwa[reads.dt[, get(col)]]

        ## in theory, all reads should align to the same strand of the contig
        ## as their sequence has been flipped so that the mates of seed reads are on the same
        ## strand as the seed read
        ## so we should count the number of positive and negative strand alignments of reads
        ## and reverse-complement contigs with more negative strand read alignments
        realns.strand.dt = as.data.table(ralns)[, .(n.positive.alns = sum(strand == "+"),
                                                    n.negative.alns = sum(strand == "-")),
                                                by = seqnames]
        rc.seqnames = realns.strand.dt[n.negative.alns > n.positive.alns,
                                       as.numeric(as.character(seqnames))]
        ## also these should be a dichotomy - there should be many more of one than the other
        good.seqnames = realns.strand.dt[pmax(n.positive.alns, n.negative.alns) >
                                         20 * pmin(n.positive.alns, n.negative.alns),
                                         as.numeric(as.character(seqnames))]
        all.tigs[rc.seqnames] = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(all.tigs[rc.seqnames])))
        all.tigs = all.tigs[good.seqnames]
            
        ## if there are contigs where reads align to the negative strand of that contig
        ## get the reverse complement
        unaln.reads = which(!(as.character(1:reads.dt[,.N]) %in% mcols(ralns)[, "qname"]))
        ## if (verbose) { message("Number of unaligned reads: ", length(unaln.reads)) }
        if (length(unaln.reads) > 3) {
            unaln.qnames = reads.dt[unaln.reads, qname]
            unaln.reads = which(reads.dt[, qname] %in% unaln.qnames)
            unaln.tigs = RSeqLib::Fermi(reads = reads.dt[unaln.reads, get(col)],
                                        qual = reads.dt[unaln.reads, get(qcol)], assemble = TRUE)
            ## if (verbose) { message("Number of new contigs: ", length(RSeqLib::contigs(unaln.tigs))) }
            if (length(RSeqLib::contigs(unaln.tigs))) {
                all.tigs = c(all.tigs, RSeqLib::contigs(unaln.tigs))
            }
        }
    }
    return(all.tigs)
}

#' @name qc_single_contig
#' @title qc_single_contig
#'
#' @description
#'
#' Takes a single contig (one qname only) and the seed locus from which the contig was assembled.
#' Annotates the sequences to which the contig was aligned (e.g. human/viral/etc.)
#' Identify whether the contig represents an ALT allele
#'
#' Potential ALT alleles include:
#' - contigs not aligning to the seed region at all
#' - contigs aligning to the seed region but with a substantial number of unmapped bases
#' - contigs representing junctions, including INV, TRA, DEL, DEUP
#' - contigs representing phased complex rearrangements
#'
#' @param calns (data.table) containing a single contig and its alignments. This must be unique.
#' @param seed.gr (GRanges) genomic range of seed region for contig assembly
#' @param low.mappability.gr (GRanges) GRanges of low mappability bases
#' @param unassembled.gr (GRanges) GRanges of unassembled bases
#' @param athresh (numeric) minimum number of aligned bases a valid alignment (default 20 bp)
#' @param seed.pad (numeric) a distal alignment must be at least this many bases away from the seed. should be close to maximum plausible insert size (default 1000 bp)
#' @param maqp.thresh (numeric) default 60, minimum alignment MAPQ to be considered a 'high MAPQ' alignment
#'
#' @return data.table with the following added columns
#' - outside.seed (logical)
#' - outside.stranded.seed (logical)
#' - alength (numeric) : total length of alignment
#' - single.chunk (logical)
#' - keep (logical) - this indicates whether the contig should be kept for further analysis
#' - fbi (logical) - is the contig a foldback-inversion
#' - unmapped bases (logical) does the contig have > athresh unmapped bases?
#' - junction (logical) does the contig represent a junction?
#' - complex (logical) does the contig represent a phased complex rearrangement?
#' - homology (numeric) number of base pairs of breakend homology (only applicable to junctions)
#' - insertion (numeric) number of inserted base pairs between breakends (only applicable to junctions)
#' - high.mapq (logical)
qc_single_contig = function(calns.dt,
                            seed.gr,
                            low.mappability.gr = GRanges(),
                            unassembled.gr = GRanges(),
                            athresh = 20,
                            seed.pad = 1e3,
                            mapq.thresh = 60)
{

    qn = calns.dt[, unique(qname)]
    if (length(qn) != 1) {
        stop("Only a single unique qname can be supplied")
    }

    
    
    ## what parts of the contig are completely outside of the seed region?
    calns.dt[, outside.seed := NA]
    calns.dt[, outside.stranded.seed := NA]
    if (calns.dt[!is.na(seqnames) & end >= start, .N]) {
        calns = dt2gr(calns.dt[!is.na(seqnames) & end >= start, .(seqnames, start, end, strand)])
        calns.dt[, outside.seed := !(calns %^% (seed.gr + seed.pad))]
        calns.dt[, outside.stranded.seed := !(calns %^^% (seed.gr + seed.pad))]
    }

    ## use countCigar to get the number of fully matching bases
    ## if alignment is valid, then count number of fully matching bases
    calns.dt[, alength := 0]
    if (calns.dt[!is.na(cigar), .N]) {
        calns.dt[!is.na(cigar), alength := bamUtils::countCigar(cigar)[, "M"]]
    }

    ## check if is single chunk and outside seed
    calns.dt[, single.chunk := all(abs(qwidth - alength) < athresh), by = qname]

    ## check if it is unaligned
    calns.dt[, unaligned := is.na(seqnames) | end < start | is.na(cigar)]

    ## decide whether to keep single chunk contigs based on whether they align outside the seed region
    keep = FALSE
    if (calns.dt[(single.chunk), .N]) {
        keep = any(calns.dt[, outside.seed], na.rm = TRUE)
    }

    ## otherwise, if the alignment is broken into multiple chunks
    ## there are three cases:
    ## - all parts of the contig have alignments within the seed region, and on the same strand --> discard
    ## - all parts of the contig have alignments within the seed region, and but on different strands --> fbi
    ## - some parts of the contig align outside the seed region --> chimeric

    ## define defaults for added annotations
    
    fbi = FALSE
    unmapped.bases = FALSE
    junction = FALSE
    jstring = NA_character_
    proximal.breakend = NA_character_
    distal.breakend = NA_character_
    complex = NA
    homology = NA
    insertion = NA
    high.mapq = NA
    decoy = FALSE
    polya = FALSE
    unassembled = FALSE
    viral = FALSE
    repetitive = FALSE
    human = FALSE
    proximal.unmappable = NA
    distal.unmappable = NA
    proximal.unassembled = NA
    distal.unassembled = NA
    distal.foreign = NA ## does the distal part of the contig align to foreign sequence?

    ## decide whether the proximal breakend is unmappable
    proximal.unmappable = seed.gr %o% low.mappability.gr

    ## decide whether proximal breakend is unassembled
    proximal.unassembled = seed.gr %o% unassembled.gr

    ## classify single chunk contigs with alignments outside the seed region
    if (calns.dt[(single.chunk) & !(unaligned), .N] & keep) {

        distal.gr = dt2gr(calns.dt[, .(seqnames, start, end, strand)])
        
        if (calns.dt[(single.chunk) & !(unaligned), .N] == 1) {
            junction = TRUE
            proximal.breakend.gr = gr.flipstrand(gr.end(seed.gr))
            proximal.breakend = paste0(seqnames(proximal.breakend.gr),
                                       ":",
                                       GenomicRanges::end(proximal.breakend.gr),
                                       strand(proximal.breakend.gr))
            distal.breakend.gr = gr.start(distal.gr)[1]
            distal.breakend = paste0(seqnames(distal.breakend.gr),
                                     ":",
                                     GenomicRanges::start(distal.breakend.gr),
                                     strand(distal.breakend.gr))
            jstring = paste0(proximal.breakend, ",", distal.breakend)
        }

        ## decide whether distal breakend is unmappable
        distal.unmappable = sum(distal.gr %o% low.mappability.gr,
                                na.rm = TRUE)

        ## decide whether distal breakend is unassembled
        distal.unassembled = sum(distal.gr %o% unassembled.gr,
                                na.rm = TRUE)

        distal.foreign = any(!(as.character(seqnames(distal.gr)) %in% as.character(seqnames(low.mappability.gr))), na.rm = TRUE)
    }

    if (calns.dt[!(single.chunk) & !(unaligned), .N]) {
        qwidth = calns.dt[qname == qn][, qwidth][1]
        ## generate gChain from contig cigar string
        gr = GRanges(seqnames = calns.dt[qname == qn, seqnames],
                     ranges = IRanges(start = calns.dt[qname == qn, start],
                                      end = calns.dt[qname == qn, end]),
                     strand = calns.dt[qname == qn, strand],
                     cigar = calns.dt[qname == qn, cigar],
                     qname = qn)
        cg.contig = gChain::cgChain(cigar = gr)
        
        ## identify which parts of the contig fail to overlap the stranded seed
        y = cg.contig$y
        x = cg.contig$x

        ## reorder each by start site on the contig
        x.order = order(GenomicRanges::start(x), decreasing = FALSE)
        x = x[x.order]
        y = y[x.order]

        ## identify check whether > athresh bases fail to align anywhere
        unmapped.bases.indicator = sum(width(x)) + athresh < qwidth

        
        ## check alignment to *STRAND-SPECIFIC* seed
        outside.stranded.seed.indicator = calns.dt[qname == qn, any(outside.stranded.seed & alength > athresh)]
        if (calns.dt[qname == qn & (!outside.stranded.seed), .N]) {
            mcols(y)[, "outside.stranded.seed"] = !(y %^^% (seed.gr + seed.pad))
            mcols(x)[, "outside.stranded.seed"] = mcols(y)[, "outside.stranded.seed"]
            stranded.seed.region.width = sum(width(x %Q% (!outside.stranded.seed)))
            outside.stranded.seed.indicator = qwidth > stranded.seed.region.width + athresh

            ## get the number of distal chunks outside the stranded seed
            distal.nchunks = length((y %Q% (outside.stranded.seed)))
        } else {
            distal.nchunks = length(y)
        }
        
        outside.unstranded.seed.indicator = calns.dt[qname == qn, any(outside.seed & alength > athresh)]
        if (calns.dt[qname == qn & (!outside.seed), .N]) {
            ## check alignment to *STRAND-AGNOSTIC* seed
            ## indicator is TRUE if there is an alignment outside of the *UN*stranded seed
            mcols(y)[, "outside.unstranded.seed"] = !(y %^% (seed.gr + seed.pad))
            mcols(x)[, "outside.unstranded.seed"] = mcols(y)[, "outside.unstranded.seed"]
            seed.region.width = sum(width(x %Q% (!outside.unstranded.seed)))
            outside.unstranded.seed.indicator = qwidth > seed.region.width + athresh

            ## get number of distal chunks
            distal.nchunks = length((y %Q% (outside.unstranded.seed)))
        } else {
            distal.nchunks = length(y)
        }
        
        ## indicate whether contig is noise or FBI
        fbi = (!outside.unstranded.seed.indicator) & (outside.stranded.seed.indicator) & (distal.nchunks == 1)
        keep = (fbi) | outside.unstranded.seed.indicator | unmapped.bases.indicator

        ## indicate whether contig represents a junction
        junction = (fbi) | (outside.unstranded.seed.indicator & (distal.nchunks == 1))
        
        ## indicate whether contig represents a phased complex rearrangement
        complex = outside.unstranded.seed.indicator & (distal.nchunks > 1)
        
        if (fbi | junction | complex) {
            if (fbi) {
                distal.gr = y %Q% (outside.stranded.seed)
            } else if (junction) {
                distal.gr = y %Q% (outside.unstranded.seed)
            } else if (complex) {
                distal.gr = y %Q% (outside.unstranded.seed)
            }
            ## is the distal part of the contig unmappable?
            distal.unmappable = sum(distal.gr %o% low.mappability.gr,
                                    na.rm = TRUE)

            ## decide whether distal breakend is unassembled
            distal.unassembled = sum(distal.gr %o% unassembled.gr,
                                     na.rm = TRUE)

            distal.foreign = any(!(as.character(seqnames(distal.gr)) %in% as.character(seqnames(low.mappability.gr))), na.rm = TRUE)
        }

        ## indicate whether there is an insertion
        insertion = NA
        homology = NA
        if ((junction | complex)) {
            ## insertion means there's a gap in the contig that isn't aligned to anything
            insertion = ifelse(GenomicRanges::start(x[2]) >= GenomicRanges::end(x[1]) + 1,
                               GenomicRanges::start(x[2]) - GenomicRanges::end(x[1]) - 1,
                               0)
            ## homology means there's some bases in the contig aligning to multiple spots in the genome
            homology = ifelse(GenomicRanges::end(x[1]) >= GenomicRanges::start(x[2]),
                              GenomicRanges::end(x[1]) - GenomicRanges::start(x[2]) + 1,
                              0)
            proximal.breakend = paste0(seqnames(y[1]), ":",
                                       GenomicRanges::end(y[1]),
                                       ifelse(strand(y[1]) == "+", "-", "+"))
            distal.breakend = paste0(seqnames(y[2]), ":",
                                     GenomicRanges::start(y[2]),
                                     strand(y[2]))
            if (junction) {
                jstring = paste0(proximal.breakend, ",", distal.breakend)
            } else if (complex) {
                y.distal = y %Q% (outside.unstranded.seed)
                jstring = paste0(proximal.breakend, ",",
                                 paste(grl.string(y.distal), collapse = ","))
            }
        } 
    }

    ## is this a high mapq contig?
    high.mapq = all(calns.dt[qname == qn, mapq >= mapq.thresh], na.rm = TRUE)

    calns.dt[, ":="(keep = keep,
                    fbi = fbi,
                    unmapped.bases = unmapped.bases,
                    junction = junction,
                    jstring = jstring,
                    proximal.breakend = proximal.breakend,
                    distal.breakend = distal.breakend,
                    complex = complex,
                    homology = homology,
                    insertion = insertion,
                    high.mapq = high.mapq,
                    viral = any(c_type %like% 'viral', na.rm = TRUE),
                    polya = any(c_type %like% 'poly', na.rm = TRUE),
                    unassembled = any(c_type %like% 'unassembled', na.rm = TRUE),
                    repetitive = any(c_type %like% 'rep', na.rm = TRUE),
                    human = all(c_type %like% 'human', na.rm = TRUE),
                    decoy = any(c_type %like% 'decoy', na.rm = TRUE),
                    proximal.unmappable = proximal.unmappable,
                    proximal.unassembled = proximal.unassembled,
                    distal.unmappable = distal.unmappable,
                    distal.unassembled = distal.unassembled,
                    distal.foreign = distal.foreign)]

    return(calns.dt)
}

#' @name qc_contigs
#' @title qc_contigs
#'
#' @description
#'
#' Given (possibly multiple) contigs and their alignments
#' Identifies contigs corresponding with ALT alleles and removes
#'
#' @param calns (data.table) alignment of possibly multiple contigs
#' @param seed.gr (GRanges) genomic ranges of seed region used to assembled contigs
#' @param low.mappability.gr (GRanges) ranges of low mappability bases
#' @param unassembled.gr (GRanges) ranges of unassembled bases
#' @param athresh (numeric) # matching bases needed for alignment
#' @param seed.pad (numeric) need at least this many bps from seed region to be considered a "distal" aln
#'
#' @return data.table.
qc_contigs = function(calns,
                      seed.gr,
                      low.mappability.gr = GRanges(),
                      unassembled.gr = GRanges(),
                      athresh = 20,
                      seed.pad = 1e3)
{
    calns.dt = copy(calns)
    if (!calns.dt[, .N]) {
        return(calns.dt)
    }
    qns = unique(calns.dt[, qname])
    ## create data table of classified contigs
    clf = lapply(qns, function(qn) {
        return(qc_single_contig(calns.dt[qname == qn,],
                                seed.gr,
                                athresh = athresh,
                                seed.pad = seed.pad,
                                low.mappability.gr = low.mappability.gr,
                                unassembled.gr = unassembled.gr))
    })
    clf = rbindlist(clf)
    return(clf)

}

#' @name find_telomeres
#' @title find_telomeres
#'
#' @param seq (character) vector of sequences
#' @param gorc (either g, c, or both, input to eighteenmer)
#' @param verbose (logical) default FALSE
#'
#' @return logical with length of seq, indicating whether that sequence contains a c or g telomere
find_telomeres = function(seq = character(), gorc = "g", verbose = FALSE) {

    if (!gorc %in% c("g", "c", "both")) {
        stop("Invalid value: gorc must be one of c('g', 'c', 'both')")
    }

    if (!length(seq)) {
        return(logical())
    }

    if (verbose) { message("searching for ", gorc, " telomeres in ", length(seq), " sequences") }

    telomere.query = eighteenmer(gorc = gorc)
    telomere.subject = Biostrings::DNAStringSet(x = seq)
    res = vwhichPDict(pdict = telomere.query, subject = telomere.subject, min.mismatch = 0, max.mismatch = 0)
    return(base::lengths(res) > 0)
}
