#' @importFrom stringr str_count

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

#' @name label_contigs
#' @title label_contigs
#'
#' @description
#'
#' Labels the category of each contig associated with a loose end, as well as the overall
#' category of the loose end. (e.g. T0/T1/T2/MYS)
#' 
#' @param li (data.table) data.table coercible to GRanges representing loose end
#' @param bps (data.table) data.table representing breakpoints coercible to GRanges. bps needs to have columns tumor.count and normal.count representing supporting barcodes for that contig.
#' @param mapq.gr (GRanges) tiled granges with field mapq0.frac
#' @param unmappable.gr (GRanges) ranges that are considered CN unmappable
#' @param mappability.thresh (numeric) 0.1 for mappable categorization
#' @param specificity.thresh (numeric) cutoff for tumor-specific contigs
#'
#' @param returns a list with two items:
#' - call (data.table with a single row and metadata column call)
#' - breakpoints (data.table labeled breakpoints. empty table if no contigs)
label_contigs = function(le, bps, mapq.gr, unmappable.gr,
                         mappability.thresh = 0.5,
                         specificity.thresh = 0) {

    ## identify the proximal category of the loose end
    ## in addition, consider the mappability of the actual loose end
    le.strongly.unmappable = dt2gr(le) %^% unmappable.gr
    if (dt2gr(le) %^% mapq.gr) {
        le.mappable = (dt2gr(le) %$% mapq.gr)$mapq0.frac <= mappability.thresh
    } else {
        le.mappable = FALSE
    }

    pcat = ifelse(le.mappable, "M", ifelse(le.strongly.unmappable, "S", "W"))

    if (!bps[, .N]) {
        call = copy(le)[, ":="(category = "MYS", proximal = pcat, distal = NA, nonspecific = NA,
                               c_telomere = FALSE, g_telomere = FALSE)]
        return(list(call = call, breakpoints = bps))
    }

    ## check whether the breakend contains contigs with telomeric sequence
    c.telo = FALSE
    if ("query_c_telomere" %in% names(bps)) {
        c.telo = any(bps[, query_c_telomere], na.rm = TRUE)
    }
    g.telo = FALSE
    if ("query_g_telomere" %in% names(bps)) {
        g.telo = any(bps[, query_g_telomere], na.rm = TRUE)
    }
    
    ## overlap candidate breakpoints with mappability masks
    bps.gr = dt2gr(bps) %$% mapq.gr
    mcols(bps.gr)[, "cn.unmappable"] = bps.gr %^% unmappable.gr
    bps = as.data.table(bps.gr)

    ## get rid of any not tumor-specific contigs
    bps[, specific := (normal.count <= specificity.thresh * tumor.count) & (tumor.count > 1)]

    ## if there are zero specific contigs, then category is MYS!
    if (!bps[(specific), .N]) {
        call = copy(le)[, ":="(category = "MYS", proximal = pcat, distal = NA, nonspecific = TRUE,
                               c_telomere = c.telo, g_telomere = g.telo)]
        return(list(call = call, breakpoints = bps))
    }

    ## label the distal breakpoints of each contig
    ## first if there are ANY mappable distal breakpoints, the distal breakpoint is M
    bps[(!proximal), any.mappable := any(mapq0.frac < mappability.thresh, na.rm = TRUE), by = .(name)]

    ## if ALL distant breakpoints are CN unmappable, the distal breakpoint is S
    ## then anything in between, is W
    bps[(!proximal), strongly.unmappable := all(.SD$cn.unmappable, na.rm = TRUE), by = .(name)]
    bps[(!proximal), distal.category := ifelse(any.mappable, "M", ifelse(strongly.unmappable, "S", "W"))]

    ## finally, NA everything that is not specific
    bps[(!specific), ":="(distal.category = NA)]

    ## in addition, consider the mappability of the actual loose end
    le.strongly.unmappable = dt2gr(le) %^% unmappable.gr
    if (dt2gr(le) %^% mapq.gr) {
        le.mappable = (dt2gr(le) %$% mapq.gr)$mapq0.frac <= mappability.thresh
    } else {
        le.mappable = FALSE
    }

    bps[(specific), proximal.category := ifelse(le.mappable, "M", ifelse(le.strongly.unmappable, "S", "W"))]
    bps[(!specific), ":="(proximal.category = NA)]

    bps[, proximal.summary := ifelse(any(proximal.category == "M", na.rm = TRUE), "M",
                              ifelse(all(proximal.category == "S", na.rm = TRUE), "S", "W"))]
    bps[, distal.summary := ifelse(any(distal.category == "M", na.rm = TRUE), "M",
                            ifelse(all(distal.category == "S", na.rm = TRUE), "S", "W"))]

    proximal = pcat##bps[, proximal.summary][1]
    distal = bps[, distal.summary][1]
    if (proximal == "M") {
        if (distal == "M") {
            call = copy(le)[, ":="(category = "T0", proximal = proximal, distal = distal, nonspecific = FALSE,
                                   c_telomere = c.telo, g_telomere = g.telo)]
        } else {
            call = copy(le)[, ":="(category = "T1", proximal = proximal, distal = distal, nonspecific = FALSE,
                                   c_telomere = c.telo, g_telomere = g.telo)]
        }
    } else {
        if (distal == "M") {
            call = copy(le)[, ":="(category = "T1", proximal = proximal, distal = distal, nonspecific = FALSE,
                                   c_telomere = c.telo, g_telomere = g.telo)]
        } else {
            call = copy(le)[, ":="(category = "T2", proximal = proximal, distal = distal, nonspecific = FALSE,
                                   c_telomere = c.telo, g_telomere = g.telo)]
        } 
    }

    return(list(call = call, breakpoints = bps))
}

#' @name check_contig_support_wrapper
#' @title check_contig_support_wrapper
#'
#' @description
#'
#' Given a set of filtered contigs, checks support for each contig using the appropriate function
#' e.g. if the contig is chimeric, uses check_split_contig_support
#' whereas, if the contig contains the distal side only, uses check_distal_only_contig_support
#'
#' Produces a data table with the name of each contig and the number of tumor/normal supporting qnames
#'
#' @param calns (data.table) alignment for POSSIBLY MULTIPLE contigs
#' @param reads.dt (data.table)
#' @param ref (BWA)
#' @param name.field (character) default "name"
#' @param verbose (logical) default FALSE
check_contig_support_wrapper = function(calns, reads.dt, ref, name.field = "name", verbose = FALSE) {

    ## get unique qnames
    if (!name.field %in% names(calns)) {
        stop("Invalid name field for contigs: ", name.field)
    }
    calns = copy(calns)[, name := get(name.field)]
    qns = unique(calns[, name])

    supporting.qname.counts = lapply(qns,
                                     function(qn) {
                                         if (verbose) { message("Checking qname: ", qn) }
                                         if (all(calns[name == qn, single.chunk])) {
                                             if (verbose) { message("Single chunk alignment for distal-only") }
                                             rs = check_distal_only_contig_support(reads.dt = reads.dt,
                                                                                   calns = calns[name == qn,],
                                                                                   verbose = verbose)
                                         } else {
                                             if (verbose) {message("Split alignment for chimeric contig")}
                                             rs = check_split_contig_support(reads.dt = reads.dt,
                                                                             calns = calns[name == qn,],
                                                                             ref = ref,
                                                                             verbose = verbose)
                                         }
                                         ## count the number of tumor and normal supporting qnames
                                         if (rs[, .N]) {
                                             tumor.count = unique(rs, by = "qname")[track %like% "sample",.N]
                                             normal.count = unique(rs, by = "qname")[track %like% "control",.N]
                                         } else {
                                             tumor.count = 0
                                             normal.count = 0
                                         }
                                         if (verbose) {
                                             message("# tumor qnames: ", tumor.count)
                                             message("# normal qnames: ", normal.count)
                                         }
                                         return(data.table(name = qn,
                                                           tumor.count = tumor.count,
                                                           normal.count = normal.count))
                                     })

    return(rbindlist(supporting.qname.counts))
}


#' @name check_distal_only_contig_support
#' @title check_distal_only_contig_support
#'
#' @description
#'
#' Gets supporting reads for contigs with distal-only alignments
#' This is likely to cause false positives for split contigs so use only with distal-only contigs
#' 
#' @param calns (data.table) represents contig alignment for a SINGLE contig
#' @param reads.dt (data.table) reads corresponding to the loose end associated with this contig
#' @param ref.pad (numeric) pad (bp) around seed region to get reference alignment, default 5 kbp
#' @param seed.pad (numeric) number of base pairs around peak to search for supporting reads, default 0
#' @param verbose (logical) default FALSE
check_distal_only_contig_support = function(calns, reads.dt, ref.pad = 5e3, seed.pad = 0, verbose = FALSE) {

    if (verbose) { message("Grabbing contig peak") }
    win = parse.gr(unique(calns[, seed]))

    if (verbose) { message("Grabbing reads near peak and their mates") }
    seed.qnames = reads.dt[dt2gr(reads.dt) %^^% (win + seed.pad), qname]
    window.reads.dt = reads.dt[qname %in% seed.qnames] ## which reads correspond with that qname?
    window.reads.gr = dt2gr(window.reads.dt)

    if (verbose) { message("Grabbing reference sequence around peak") }
    refseq = Biostrings::getSeq(Hsapiens, trim(gr.fix(gr.chr(win + ref.pad), Hsapiens)))
    ref = RSeqLib::BWA(seq = as.character(refseq))

    if (verbose) { message("Building contig gChain") }
    calns.gr = dt2gr(calns)
    if (!is.null(calns.gr$query.seq)) {
        calns.gr$seq = calns.gr$query.seq
    }
    cg.contig = gChain::cgChain(calns.gr)
    
    rs = contig.support(reads = window.reads.gr,
                        contig = calns.gr,
                        ref = ref,
                        cg.contig = cg.contig,
                        verbose = verbose)

    window.reads.dt[, supporting := qname %in% rs$qname]
    return(window.reads.dt[(supporting),])
}

#' @name check_split_contig_support
#' @title check_split_contig_support
#'
#' @description
#'
#' Gets supporting reads for split contigs (aligning to two separate places in the reference)
#' Note that this is not appropriate for contigs that represent purely the mate sequence
#' 
#' Takes as input a contig alignment, reads, and a reference BWA object
#'
#' @param calns (data.table) represents contig alignment
#' @param reads.dt (data.table) reads corresponding to the loose end associated with this contig
#' @param ref (BWA) reference to compare alignment against
#' @param seed.pad (numeric) number of base pairs around peak to search for supporting reads, default 0
#' @param verbose (logical) default FALSE
check_split_contig_support = function(calns, reads.dt, ref, seed.pad = 0, verbose = FALSE) {

    if (verbose) { message("Grabbing contig peak") }
    win = parse.gr(unique(calns[, seed]))

    if (verbose) { message("Grabbing reads near peak and their mates") }
    seed.qnames = reads.dt[dt2gr(reads.dt) %^^% (win + seed.pad), qname]
    window.reads.dt = reads.dt[qname %in% seed.qnames] ## which reads correspond with that qname?
    window.reads.gr = dt2gr(window.reads.dt)

    if (verbose) { message("Building contig gChain") }
    calns.gr = dt2gr(calns)
    if (!is.null(calns.gr$query.seq)) {
        calns.gr$seq = calns.gr$query.seq
    }
    cg.contig = gChain::cgChain(calns.gr)
    
    rs = contig.support(reads = window.reads.gr,
                        contig = calns.gr,
                        ref = ref,
                        cg.contig = cg.contig,
                        verbose = verbose)

    window.reads.dt[, supporting := qname %in% rs$qname]
    return(window.reads.dt[(supporting),])
}

#' @name grab_contig_breakends_wrapper
#' @title grab_contig_breakends_wrapper
#'
#' @param calns (data.table) contig alignments
#' @param name.field (character) default name
#' @param verbose (logical)
grab_contig_breakends_wrapper = function(calns, name.field = "name", verbose = FALSE) {

    ## get unique qnames
    if (!name.field %in% names(calns)) {
        stop("Invalid name field for contigs: ", name.field)
    }
    calns = copy(calns)[, name := get(name.field)]
    qns = unique(calns[, name])

    bnds = lapply(qns, function(qn) {
        if (verbose) { message("Checking qname: ", qn) }
        seed = parse.gr(unique(calns[name == qn, seed]))
        qn.bnd = grab_contig_breakends(calns = calns[name == qn,],
                                       seed = seed)
        if (qn.bnd[, .N]) {
            qn.bnd[, name := qn]
        }
    })

    bnds = rbindlist(bnds, fill = TRUE, use.names = TRUE)

    return(bnds)
}


#' @name grab_contig_breakends
#' @title grab_contig_breakends
#'
#' @description
#'
#' Get seed side and mate side breakends for a single contig.
#'
#' Assumes that a unique qname is supplied.
#' Also assumes that the alignment is split.
#'
#' If more than one contig is supplied, throws an exception.
#' 
#' @param calns (GRanges) contig alignment
#' @param seed (GRanges) seed region of the contig
#' @param reduce.pad (numeric) reduce multipmapping regions using this pad (default 20)
#' @param seed.pad (numeric) bp pad to be considered lying outside peak (default 1e3)
#' @param verbose (logical)
#'
#' @return data.table with seqnames, start, end, strand, and a few metadata columns:
#' - proximal (logical) proximal vs. distal breakends?
#' - proximal.first (logical) proximal part of the contig is in the first part of the contig
#' - seed.only (logical) contig is seed only
#' - mate.only (logical) conig is mate only
#' - c_type (character)
#' - c_spec (character)
#' - cigar (character)
#' strand is in breakend orientation (so if there were a junction it would be consistent with the junction)
grab_contig_breakends = function(calns, seed, reduce.pad = 20, seed.pad = 1e3, verbose = FALSE) {

    empty.res = data.table(seqnames = character(), start = numeric(), end = numeric(), strand = character(),
                           proximal = logical(), proximal.first = logical(),
                           seed.only = logical(), mate.only = logical())
    
    if (!length(calns)) {
        stop("Empty GRanges supplied for contigs")
    }

    if (length(unique(calns$qname)) > 1) {
        stop("Multiple contigs supplied!")
    }

    ## extract some metadata...
    fbi = FALSE
    if (!is.null(calns$fbi)) {
        fbi = calns[, fbi][1]
        fbi = ifelse(is.na(fbi), FALSE, fbi)
    }

    ## create cgChain for this contig
    cgtigs = gChain::cgChain(calns, sn = calns$qname)

    ## extract x and y side of mappings and overlap y side with peak
    y = cgtigs$y
    mcols(y)[, "peak"] = y %^% (seed + seed.pad)
    x = cgtigs$x
    mcols(x)[, "peak"] = mcols(y)[, "peak"]

    ## if we have an inversion, then use strand-specific peak
    ## if (all(mcols(y)[, "peak"])) {
    if (fbi) {
        ## this implies that we likely have an inversion
        ## so we would like to identify which overlaps the peak in a strand specific way
        if (verbose) { message("Inversion detected? Using strand-specific peak") }
        mcols(y)[, "peak"] = y %^^% (seed + seed.pad)
        mcols(x)[, "peak"] = mcols(y)[, "peak"]
        ## return(empty.res)
    }

    ## if we have a mate only contig
    if (!any(mcols(y)[, "peak"])) {
        y.red = reduce(y + reduce.pad) - reduce.pad
        res = as.data.table(gr.start(y.red))[, ":="(proximal = FALSE, proximal.first = NA,
                                                    seed.only = FALSE, mate.only = TRUE)]
    } else if (all(mcols(y)[, "peak"])) {
        y.red = reduce(y + reduce.pad) - reduce.pad
        res = as.data.table(gr.start(y.red))[, ":="(proximal = TRUE, proximal.first = NA,
                                                    seed.only = TRUE, mate.only = FALSE)]
    } else {

        ## reduce the peak and non-peak regions of the alignment
        ## after adding a pad
        y.peak = reduce((y %Q% (peak)) + reduce.pad) - reduce.pad
        y.nonpeak = reduce((y %Q% (!peak)) + reduce.pad) - reduce.pad

        ## check whether the peak side is the front or back of the contig
        x.red = reduce((x %Q% (peak)) + reduce.pad) - reduce.pad

        if (start(x.red) < 20) {
            ## if the peak side is at the front of the contig
            ## the proximal bp is at the end of the peak region
            ## and the strand must be flipped
            proximal.bp = gr.flipstrand(gr.end(y.peak))
            ## the distal bp is at the start of the non peak region
            distal.bp = gr.start(y.nonpeak)
            proximal.start = TRUE
        } else {
            ## otherwise proximal is the start of the peak region
            proximal.bp = gr.start(y.peak)
            ## and distal is the end of the non-peak region with strand flipped
            distal.bp = gr.flipstrand(gr.end(y.nonpeak))
            proximal.start = FALSE
        }

        res = rbind(as.data.table(proximal.bp)[, ":="(proximal = TRUE,
                                                      proximal.first = proximal.start,
                                                      seed.only = FALSE,
                                                      mate.only = FALSE)],
                    as.data.table(distal.bp)[, ":="(proximal = FALSE,
                                                    proximal.first = proximal.start,
                                                    seed.only = FALSE,
                                                    mate.only = FALSE)])
    }

    ## add back some postentially useful metadata
    if (res[, .N]) {
        extra.cols = mcols( dt2gr(res)[, c()] %$% dt2gr(calns)[, c("c_type", "c_spec", "cigar")] )
        res[, c_type := extra.cols$c_type]
        res[, c_spec := extra.cols$c_spec]
        res[, cigar := extra.cols$cigar]
        res[, seed := gr.string(seed)]
    }
    return(res)
}

#' @name build_contigs_wrapper
#' @title build_contigs_wrapper
#'
#' @param gr (GRanges) which breakpoint do we care about
#' @param reads.dt (data.table)
#' @param ref (BWA object)
#' @param window (numeric) how many bases around window to consider for assembly (default 5e3)
#' @param assembly.region (numeric) how many bases at a time to consider for assembly (default 1e3)
#' @param stride (numeric) stride for tiling window (default 500)
build_contigs_wrapper = function(gr, reads.dt, ref,
                                 window = 5e3,
                                 assembly.region = 1e3,
                                 stride = 500,
                                 verbose = FALSE) {
    tiles = unlist(GenomicRanges::slidingWindows(x = gr + window,
                                                 width = assembly.region,
                                                 step = stride))

    all.contigs = lapply(seq_along(tiles),
                  function(ix) {
                      if (verbose) {message("Starting analysis for window: ", ix, " of ", length(tiles))}
                      seed.frame.dt = grab_seed_frame(reads.dt,
                                                      seed.gr = tiles[ix],
                                                      seq.field = "reading.frame")
                      if (seed.frame.dt[, .N]) {
                          ctigs = build_contigs(seed.frame.dt, verbose = verbose)
                          if (length(ctigs)) {
                              if (verbose) { message("Aligning contigs to reference") }
                              aln.ctigs = align_contigs(ctigs, ref, verbose = verbose, keep.unaligned = FALSE)
                              ## aln.ctigs = ref[ctigs]
                              if (verbose) { message("Filtering contigs based on structure") }
                              qc.ctigs = qc_contigs(aln.ctigs, tiles[ix])
                              if (qc.ctigs[, .N]) {
                                  qc.ctigs[, name := paste(ix, qname, sep = ".")]
                                  qc.ctigs[, seed := gr.string(tiles[ix])]
                                  ## add back the original sequence!
                                  ## qc.ctigs[, query.seq := ctigs[as.numeric(qname)]]
                              }
                          } else {
                              qc.ctigs = data.table()
                          }
                      } else {
                          qc.ctigs = data.table()
                      }
                      return(qc.ctigs)
                  })

    all.contigs = rbindlist(all.contigs, fill = TRUE, use.names = TRUE)
    return(all.contigs)
}

#' @name align_contigs
#' @title align_contigs
#'
#' @descripion
#'
#' Align contigs (character vector) to reference
#' Then aannotate the segments of the alignment according to reference seqnames
#'
#' Relies on seqnames --> annotation mapping included as part of the package...
#'
#' @param tigs (character)
#' @param ref (RSeqLib::BWA)
#' @param refseq.fn (character) path to seqnames annotation file (included)
#' @param primary.only (logical) default TRUE
#' @param remove.decoy (logical) exclude alignment to decoy sequences? default TRUE
#' @param keep.unaligned (logical) keep unaligend contigs?? default FALSE
#' @param verbose (logical) default FALSE
#'
#' @return data.table with contig alignment and metadata fields c_type and c_spec
align_contigs = function(tigs, ref,
                         refseq.fn = system.file("extdata", "reference.sequences.rds", package = "loosends"),
                         primary.only = TRUE,
                         remove.decoy = TRUE,
                         keep.unaligned = FALSE,
                         verbose = FALSE) {

    refseq.dt = readRDS(refseq.fn)

    if (verbose) { message("Finding reference alignments") }
    
    aln.tigs = ref[tigs]
    
    if (length(aln.tigs)) {
        mcols(aln.tigs)[, "c_type"] = refseq.dt[as.character(seqnames(aln.tigs)), c_type]
        mcols(aln.tigs)[, "c_spec"] = refseq.dt[as.character(seqnames(aln.tigs)), c_spec]
        mcols(aln.tigs)[, "query.seq"] = tigs[as.numeric(mcols(aln.tigs)[, "qname"])]

        if (verbose) { message("Searching alignments for C and G telomeres") }
        ## ## does the aligned chunk specifically have a C/G telomere?
        mcols(aln.tigs)[, "aln_c_telomere"] = find_telomeres(seq = mcols(aln.tigs)[, "seq"], gorc = "c")
        mcols(aln.tigs)[, "aln_g_telomere"] = find_telomeres(seq = mcols(aln.tigs)[, "seq"], gorc = "g")
        ## ## does the query overall have a C/G telomere?
        ## mcols(aln.tigs)[, "query_c_telomere"] = find_telomeres(seq = mcols(aln.tigs)[, "query.seq"], gorc = "c")
        ## mcols(aln.tigs)[, "query_g_telomere"] = find_telomeres(seq = mcols(aln.tigs)[, "query.seq"], gorc = "g")
        
        
    }

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

    if (aln.tigs.dt[, .N]) {

        if (verbose) { message("Searching query sequences for C and G telomeres") }
        
        aln.tigs.dt[, query_c_telomere := find_telomeres(seq = mcols(aln.tigs)[, "query.seq"], gorc = "c")]
        aln.tigs.dt[, query_g_telomere := find_telomeres(seq = mcols(aln.tigs)[, "query.seq"], gorc = "g")]
    }

    if (aln.tigs.dt[, .N]) {
        if (verbose) { message("Annotating primary alignments") }
        aln.tigs.dt[, primary := bamUtils::bamflag(flag)[, "isNotPrimaryRead"] == 0]
        if (primary.only) {
            aln.tigs.dt = aln.tigs.dt[(primary),]
        }
        if (verbose) {
            message("Checking for alignments to decoy sequences")
        }
        aln.tigs.dt[, decoy := as.character(seqnames) %in% refseq.dt[c_type == "decoy", seqnames]]
        if (remove.decoy) {
            aln.tigs.dt = aln.tigs.dt[(!decoy)]
        }
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
#' @param max.n (numeric) maximum number of low-quality "N" bases (default 0)
#' @param verbose (logical) default FALSE
grab_seed_frame = function(reads.dt, seed.gr, seq.field = "reading.frame", max.n = 0, verbose = FALSE) {

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

    ## for non-seed reads, reverse complement
    valid.reads.dt[, rc.reading.frame := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(get(seq.field))))]
    valid.reads.dt[, seed.frame := ifelse(seed, get(seq.field), rc.reading.frame)]
    ## valid.reads.dt[, seed.frame := ifelse(seed , get(seq.field), rc.reading.frame)]
    ## valid.reads.dt[, seed.frame := ifelse((seed & stranded.seed) | ((!seed) & (!stranded.seed)) ,
    ##                                       get(seq.field),
    ##                                       rc.reading.frame)]

    ## get rid of anything with too many N's or zero length strings
    valid.reads.dt = valid.reads.dt[stringr::str_count(get(seq.field), "N") <= max.n & nchar(get(seq.field)) > 0,]

    return(valid.reads.dt)
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

    if (!reads.dt[, .N]) {
        return(character())
    }

    ## TRY THREE CONDITIONS
    ## assembly everything
    if (verbose) {
        message("Trying assembly with all reads")
    }
    tigs = RSeqLib::Fermi(reads = reads.dt[, get(col)], qual = reads.dt[, get(qcol)], assemble = TRUE)
    all.tigs = RSeqLib::contigs(tigs)

    ## assemble discordant pairs only
    if (reads.dt[(!concord), .N] > 4) {
        if (verbose) {message("Trying assembly with only discordant reads")}
        discordant.tigs = RSeqLib::Fermi(reads = reads.dt[(!concord), get(col)],
                                         qual = reads.dt[(!concord), get(qcol)],
                                         assemble = TRUE)
        all.tigs = c(all.tigs, RSeqLib::contigs(discordant.tigs))
    }

    ## assemble loose pairs only
    if (reads.dt[(loose.pair), .N] > 4) {
        ## browser()
        if (verbose) {message("Trying assembly with only loose reads")}
        loose.tigs = RSeqLib::Fermi(reads = reads.dt[(loose.pair), get(col)],
                                    qual = reads.dt[(loose.pair), get(qcol)],
                                    assemble = TRUE)
        all.tigs = c(all.tigs, RSeqLib::contigs(loose.tigs))
    }

    ## then align reads to all of the contigs found so far
    ## if there are unaligned reads, retry assembly just from the unaligned reads
    if (length(all.tigs)) {
        if (verbose) { message("Checking for reads that don't align to any contigs") }
        tigs.bwa = BWA(seq = all.tigs)
        ralns = tigs.bwa[reads.dt[, get(col)]]
        unaln.reads = which(!(as.character(1:reads.dt[,.N]) %in% mcols(ralns)[, "qname"]))
        if (verbose) { message("Number of unaligned reads: ", length(unaln.reads)) }
        if (length(unaln.reads) > 3) {
            unaln.qnames = reads.dt[unaln.reads, qname]
            unaln.reads = which(reads.dt[, qname] %in% unaln.qnames)
            unaln.tigs = RSeqLib::Fermi(reads = reads.dt[unaln.reads, get(col)],
                                        qual = reads.dt[unaln.reads, get(qcol)], assemble = TRUE)
            if (verbose) { message("Number of new contigs: ", length(RSeqLib::contigs(unaln.tigs))) }
            if (length(RSeqLib::contigs(unaln.tigs))) {
                all.tigs = c(all.tigs, RSeqLib::contigs(unaln.tigs))
            }
        }
    }
    return(all.tigs)
}

#' @name qc_contigs
#' @title qc_contigs
#'
#' @description
#'
#' Given contigs and their alignments, remove contigs that are uninformative for loose end classification
#'
#' @param calns (data.table) alignment of possibly multiple contigs
#' @param athresh (numeric) # matching bases needed for alignment
#' @param seed.pad (numeric) need at least this many bps from seed region to be considered a "distal" aln
#'
#' @return data.table.
#' keeps any metadata that was already a part of calns, and adds the following columns
#' - outside.seed (logical)
#' - outside.stranded.seed (logical)
#' - alength (numeric) : total length of alignment
#' - single.chunk (logical)
#' - keep (logical) - this indicates whether the contig should be kept for further analysis
#' - fbi (logical) - is the contig a foldback-inversion
#' - unmapped bases (logical) does the contig have > athresh unmapped bases?
qc_contigs = function(calns, seed.gr, athresh = 20, seed.pad = 1e3) {

    calns.dt = copy(calns)

    if (!calns.dt[, .N]) {
        return(calns.dt)
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

    ## keep unaligned contigs if they are big
    if (calns.dt[(unaligned), .N]) {
        calns.dt[(unaligned), keep := qwidth > athresh]
    }

    ## for alignments that are a single chunk,
    ## all such alignments need to be outside the seed region
    ## or all alignments need to have at least 20 bps that do not match
    if (calns.dt[(single.chunk) & (!unaligned), .N]) {
        calns.dt[(single.chunk), keep := all(outside.seed) & any(alength > athresh), by = qname]
        ## add other metadata columns which usually pertain only to split contigs
        calns.dt[(single.chunk), fbi := FALSE]
        calns.dt[(single.chunk), unmapped.bases := FALSE]
    }

    ## otherwise, if the alignment is broken into multiple chunks
    ## there are three cases:
    ## - all parts of the contig have alignments within the seed region, and on the same strand --> discard
    ## - all parts of the contig have alignments within the seed region, and but on different strands --> fbi
    ## - some parts of the contig align outside the seed region --> chimeric

    ## we can distinguish between these three cases using cgChain
    if (calns.dt[!(single.chunk) & !(unaligned), .N]) {

        ## iterate over unique qnames
        qns = unique(calns.dt[(!single.chunk) & (!unaligned), qname])

        clf = lapply(qns,
                     function (qn) {
                         ## get width of the contig
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
                         ## identify check whether > athresh bases fail to align anywhere
                         unmapped.bases.indicator = sum(width(x)) + athresh < qwidth
                         ## check alignment to *STRAND-SPECIFIC* seed
                         outside.stranded.seed.indicator = calns.dt[qname == qn, any(outside.stranded.seed & alength > athresh)]
                         if (calns.dt[qname == qn & (!outside.stranded.seed), .N]) {
                             mcols(y)[, "outside.stranded.seed"] = !(y %^^% (seed.gr + seed.pad))
                             mcols(x)[, "outside.stranded.seed"] = mcols(y)[, "outside.stranded.seed"]
                             stranded.seed.region.width = sum(width(x %Q% (!outside.stranded.seed)))
                             outside.stranded.seed.indicator = qwidth > stranded.seed.region.width + athresh
                         }
                         outside.unstranded.seed.indicator = calns.dt[qname == qn, any(outside.seed & alength > athresh)]
                         if (calns.dt[qname == qn & (!outside.seed), .N]) {
                             ## check alignment to *STRAND-AGNOSTIC* seed
                             ## indicator is TRUE if there is an alignment outside of the *UN*stranded seed
                             mcols(y)[, "outside.unstranded.seed"] = !(y %^% (seed.gr + seed.pad))
                             mcols(x)[, "outside.unstranded.seed"] = mcols(y)[, "outside.unstranded.seed"]
                             seed.region.width = sum(width(x %Q% (!outside.unstranded.seed)))
                             outside.unstranded.seed.indicator = qwidth > seed.region.width + athresh
                         }
                         ## indicate whether contig is noise or FBI
                         fbi = (!outside.unstranded.seed.indicator) & (outside.stranded.seed.indicator)
                         keep = (fbi) | outside.unstranded.seed.indicator | unmapped.bases.indicator
                         res = data.table(qname = qn,
                                          fbi = fbi,
                                          unmapped.bases = unmapped.bases.indicator,
                                          keep = keep)
                         return(res)
                     })


        clf = rbindlist(clf)
        setkey(clf, "qname")
        
        calns.dt[qname %in% qns, keep := clf[qname, keep]]
        calns.dt[qname %in% qns, fbi := clf[qname, fbi]]
        calns.dt[qname %in% qns, unmapped.bases := clf[qname, unmapped.bases]]                 
     }
    return(calns.dt)
}
