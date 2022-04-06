#' @name call_loose_end2
#' @title call_loose_end2
#'
#' @description
#'
#' labels a single loose end
#'
#' @param li (data.table, GRanges, or character coercible to GRanges)
#' @param ri (reads from loosereads_wrapper)
#' @param id (character) sample id
#' @param concat.bwa (RSeqLib::BWA) reference BWA object with foreign sequences
#' @param human.bwa (RSeqLib::BWA) reference BWA with human only sequences
#' @param window (numeric) window in BP for building contigs, default 5e3
#' @param low.mappability.gr.fn (character) path to low mappability ranges
#' @param unassembled.gr.fn (character) path to low mappability ranges
#' @param verbose (logical) default FALSE
#' @param ... (additional inputs to various helper functions)
#'
#' @return list with four items
#' - call
#' - contigs (curated contigs)
#' - distal.bps (distal breakends)
#' - all.bps (all breakends)
#' - all.tigs (all contigs)
call_loose_end2 = function(li,
                           ri,
                           id,
                           concat.bwa,
                           human.bwa,
                           window = 5e3,
                           low.mappability.gr.fn = system.file("extdata", "wide.multimapping.ranges.rds",
                                                               package = "loosends"),
                           unassembled.gr.fn = system.file("extdata", "assembly.gaps.ranges.rds",
                                                           package = "loosends"),
                           verbose = FALSE,
                           ...) {

    args = list(...)

    if (verbose) {message("Preparing loose reads")}
    le.dt = prep_loose_ends(li = li, id = id)
    prepped.reads.dt = prep_loose_reads(li = le.dt, loose.reads.dt = ri)

    if (verbose) { message("Reading multi-mapping ranges") }
    ## read unmappable sequences
    if (file.exists(low.mappability.gr.fn) && file.info(low.mappability.gr.fn)$size > 0) {
        low.mappability.gr = readRDS(low.mappability.gr.fn)
        if (!inherits(low.mappability.gr, 'GRanges')) {
            stop("low.mappability.gr.fn must contain GRanges")
        }
    }

    ## read unassembled ranges
    if (file.exists(unassembled.gr.fn) && file.info(unassembled.gr.fn)$size > 0) {
        unassembled.gr = readRDS(unassembled.gr.fn)
        if (!inherits(unassembled.gr, 'GRanges')) {
            stop("unassembled.gr.fn must contain GRanges")
        }
    }

    if (verbose) {message("Building contigs")}
    all.tigs = build_contigs_wrapper(gr = gr.flipstrand(dt2gr(le.dt)),
                                     reads.dt = prepped.reads.dt,
                                     ref = concat.bwa,
                                     window = window,
                                     low.mappability.gr = low.mappability.gr,
                                     unassembled.gr = unassembled.gr,
                                     verbose = verbose)

    keep.tigs = copy(all.tigs)
    if (all.tigs[, .N]) {
        keep.tigs = all.tigs[(keep),]
    }

    if (verbose) { message("Number of potential ALT contigs based on structure: ", keep.tigs[, .N]) }

    keep.tigs.support = check_contig_support_wrapper(reads.dt = prepped.reads.dt,
                                                     calns = keep.tigs,
                                                     ref = human.bwa,
                                                     verbose = verbose)

    label.dt = check_contig_concordance(keep.tigs.support, verbose = verbose)

    ## add back loose end to help with analysis/merging
    label.dt[, loose.end := paste0(le.dt[, seqnames], ":", le.dt[, start], le.dt[, strand])]

    if (keep.tigs.support[, .N]) {
        keep.tigs.support[, loose.end := paste0(le.dt[, seqnames],
                                                ":",
                                                le.dt[, start],
                                                le.dt[, strand])]
    }

    ## make one unique table
    return(list(all.tigs = all.tigs,
                keep.tigs.support = keep.tigs.support,
                label.dt = label.dt))
    
}

#' @name check_contig_concordance
#' @title check_contig_concordance
#'
#' @param calns (data.table of contig alignments) needs to have columns "tumor.specific" and "tumor.support"
#' @param verbose (logical) default FALSE
#'
#' @return data.table with one row
#' with the following columns
#' - any.alt (logical)
#' - somatic.alt (logical)
#' - junction (logical) are any of the contigs high-mapq junctions?
#' - complex (logical) are any of the contigs high-mapq phased complex rearrangements (no junctions)?
#' - single breakend (logical) no high-mapq junctions or phased complex rearrangements, but there are ALT contigs.
#' - insertion (logical) do the junctions contain insertions? (NA if not a junction)
#' - homology (logical) do the junctions contain nonzero breakend homology? (NA if not a junction)
#' - telomeric (logical) do any tumor-specific contigs contain telomeric sequence?
#' - jstring (character) grl of the consensus junction (NA if not a junctions)
check_contig_concordance = function(calns, verbose = FALSE) {

    res = data.table(any.alt = FALSE,
                     somatic.alt = FALSE,
                     junction = FALSE,
                     complex = FALSE,
                     single.breakend = FALSE,
                     mystery = TRUE,
                     insertion = NA,
                     homology = NA,
                     telomeric = NA,
                     jstring = NA_character_,
                     polya = NA,
                     unassembled = NA,
                     viral = NA,
                     decoy = NA,
                     repetitive = NA,
                     human = NA)

    if (calns[, .N]) {
        if (calns[(tumor.support), .N] | calns[(normal.support), .N]) {

            res[, any.alt := TRUE]
            res[, mystery := FALSE]

            if (calns[(tumor.specific), .N]) {
                calns = calns[(tumor.specific),]
                res[, somatic.alt := TRUE]
            }

            ## defaults
            res[, junction := FALSE]
            res[, complex := FALSE]
            res[, single.breakend := FALSE]

            if (calns[, any(junction & high.mapq, na.rm = TRUE)]) {
                res[, junction := TRUE]
                res[, jstring := calns[(junction) & (high.mapq), jstring][1]]
                if (calns[(junction) & (high.mapq), any(!is.na(insertion), na.rm = TRUE)]) {
                    res[, insertion := calns[(junction) & (high.mapq), max(insertion, na.rm = TRUE)]]
                }
                if (calns[(junction) & (high.mapq), any(!is.na(homology), na.rm = TRUE)]) {
                    res[, homology := calns[(junction) & (high.mapq), max(homology, na.rm = TRUE)]]
                }
                ## add c_type
                res[, ":="(viral = any(calns[(junction) & (high.mapq), c_type %like% 'viral'], na.rm = TRUE),
                           polya = any(calns[(junction) & (high.mapq), c_type %like% 'poly'], na.rm = TRUE),
                           unassembled = any(calns[(junction) & (high.mapq),
                                                   c_type %like% 'unassembled'],
                                             na.rm = TRUE),
                           decoy = any(calns[(junction) & (high.mapq), c_type %like% 'decoy'], na.rm = TRUE),
                           human = all(calns[(junction) & (high.mapq), c_type %like% 'human'], na.rm = TRUE),
                           repetitive = any(calns[(junction) & (high.mapq), c_type %like% 'rep'], na.rm = TRUE))]

                ## add unmappable
                res[, ":="(proximal.unmappable = median(calns[(junction) & (high.mapq), proximal.unmappable],
                                                        na.rm = TRUE),
                           proximal.unassembled = median(calns[(junction) & (high.mapq), proximal.unassembled],
                                                        na.rm = TRUE),
                           distal.unmappable = median(calns[(junction) & (high.mapq), distal.unmappable],
                                                        na.rm = TRUE),
                           distal.unassembled = median(calns[(junction) & (high.mapq), distal.unassembled],
                                                        na.rm = TRUE),
                           distal.foreign = any(calns[(junction) & (high.mapq), distal.foreign],
                                                na.rm = TRUE))]
                           
                
            } else if (calns[, any(complex & high.mapq, na.rm = TRUE)]) {
                res[, complex := TRUE]
                res[, jstring := calns[(complex) & (high.mapq), jstring][1]]

                res[, ":="(viral = any(calns[(complex) & (high.mapq), c_type %like% 'viral'], na.rm = TRUE),
                           polya = any(calns[(junction) & (high.mapq), c_type %like% 'poly'], na.rm = TRUE),
                           unassembled = any(calns[(junction) & (high.mapq),
                                                   c_type %like% 'unassembled'],
                                             na.rm = TRUE),
                           decoy = any(calns[(complex) & (high.mapq), c_type %like% 'decoy'], na.rm = TRUE),
                           human = all(calns[(complex) & (high.mapq), c_type %like% 'human'], na.rm = TRUE),
                           repetitive = any(calns[(complex) & (high.mapq), c_type %like% 'rep'], na.rm = TRUE))]

                res[, ":="(proximal.unmappable = median(calns[(complex) & (high.mapq), proximal.unmappable],
                                                        na.rm = TRUE),
                           proximal.unassembled = median(calns[(complex) & (high.mapq), proximal.unassembled],
                                                        na.rm = TRUE),
                           distal.unmappable = median(calns[(complex) & (high.mapq), distal.unmappable],
                                                        na.rm = TRUE),
                           distal.unassembled = median(calns[(complex) & (high.mapq), distal.unassembled],
                                                        na.rm = TRUE),
                           distal.foreign = any(calns[(complex) & (high.mapq), distal.foreign], na.rm = TRUE))]
                
            } else {
                res[, single.breakend := TRUE]

                res[, ":="(viral = any(calns[, c_type %like% 'viral'], na.rm = TRUE),
                           decoy = any(calns[, c_type %like% 'decoy'], na.rm = TRUE),
                           polya = any(calns[, c_type %like% 'poly'], na.rm = TRUE),
                           unassembled = any(calns[, c_type %like% 'unassembled'], na.rm = TRUE),
                           human = all(calns[, c_type %like% 'human'], na.rm = TRUE),
                           repetitive = any(calns[, c_type %like% 'rep'], na.rm = TRUE))]

                res[, ":="(proximal.unmappable = median(calns[, proximal.unmappable],
                                                        na.rm = TRUE),
                           proximal.unassembled = median(calns[, proximal.unassembled],
                                                         na.rm = TRUE),
                           distal.unmappable = median(calns[, distal.unmappable],
                                                      na.rm = TRUE),
                           distal.unassembled = median(calns[, distal.unassembled],
                                                       na.rm = TRUE),
                           distal.foreign = any(calns[, distal.foreign], na.rm = TRUE))]
                
            }

            if (calns[, any(query_c_telomere, na.rm = TRUE)] |
                calns[, any(aln_c_telomere, na.rm = TRUE)] |
                calns[, any(query_g_telomere, na.rm = TRUE)] |
                calns[, any(aln_g_telomere, na.rm = TRUE)]) {

                res[, telomeric := TRUE]
                
            } else {
                res[, telomeric := FALSE]
            }
        }
    }

    ## add final annotation category
    res[, somatic := ifelse(any.alt, ifelse(somatic.alt, "somatic", "germline"), "no contigs")]
    res[, annotation := ifelse(any.alt, ifelse(junction, "junction", ifelse(complex, "complex", "breakend")), "no contigs")]

    if (verbose) {
        message("Annotation: ", res$annotation)
        message("Somatic: ", res$somatic)
    }

    return(res)
}

    ## if (verbose) { message("Grabbing contig breakends") }
    ## keep.tigs.breakends = grab_contig_breakends_wrapper(calns = keep.tigs, verbose = verbose)

    ## if (verbose) { message("Grabbing distal breakends") }
    ## distal.bps = copy(keep.tigs.breakends)
    ## final.tigs = copy(keep.tigs)
    ## if (keep.tigs.breakends[, .N]) {
    ##     setkey(keep.tigs.breakends, "name")
    ##     setkey(keep.tigs.support, "name")
    ##     keep.tigs.breakends = keep.tigs.support[keep.tigs.breakends]
    ##     distal.bps = grab_distal_breakends(keep.tigs.breakends, specificity.thresh = specificity.thresh)
    ##     final.tigs = keep.tigs[(name %in% distal.bps[, name]), ]
    ## }

    ## ## label loose end based on relation to unmappable segment
    ## proximal.mappable = !(dt2gr(le.dt) %^% (unmappable.gr + unmappable.pad))
    ## distal.mappable = NA
    ## clustered = NA
    ## if (distal.bps[, .N]) {
    ##     clustered = all(distal.bps[, clustered])
    ##     if (clustered) {
    ##         distal.mappable = any(!(dt2gr(distal.bps) %^% (unmappable.gr + unmappable.pad)))
    ##     } else {
    ##         distal.mappable = all(!(dt2gr(distal.bps) %^% (unmappable.gr + unmappable.pad)))
    ##     }
    ## }

    ## if (verbose) {
    ##     message("Proximal contig mappable: ", proximal.mappable)
    ##     message("Distal contig mappable: ", distal.mappable)
    ##     message("Distal contig breakends clustered: ", clustered)
    ## }

    ## ## decide whether any contigs align to telomeres
    ## aln.g.telomere = FALSE
    ## aln.c.telomere = FALSE
    ## if (final.tigs[, .N]) {
    ##     aln.g.telomere = final.tigs[, any(aln_g_telomere)]
    ##     aln.c.telomere = final.tigs[, any(aln_c_telomere)]
    ## }

    ## if (verbose) {
    ##     message("Contig contains G telomere: ", aln.g.telomere)
    ##     message("Contig contains C telomere: ", aln.c.telomere)
    ## }

    ## ## develop final table
    ## category = ifelse(proximal.mappable,
    ##                   ifelse(distal.mappable, "T0", "T1"),
    ##                   ifelse(distal.mappable, "T1", "T2"))
    ## if (is.na(distal.mappable)) {
    ##     category = "MYS"
    ## }

    ## if (verbose) { message("Category: ", category) }

    ## call.dt = data.table(loose.end = paste0(le.dt[, seqnames], ":", le.dt[, start], le.dt[, strand]),
    ##                      alt.contigs = keep.tigs[, .N] > 0, ## any ALT?
    ##                      supp.contigs = ifelse(keep.tigs.breakends[, .N] > 0,
    ##                                            keep.tigs.breakends[, any(tumor.count > 0 | normal.count> 0)],
    ##                                            FALSE), ## any contigs with ANY support?
    ##                      tumor.spec = final.tigs[, .N] > 0, ## any tumor-specific ALT?
    ##                      sample = id,
    ##                      proximal.mappable = proximal.mappable,
    ##                      distal.mappable = distal.mappable,
    ##                      clustered = clustered,
    ##                      aln.g.telomere = aln.g.telomere,
    ##                      aln.c.telomere = aln.c.telomere,
    ##                      category = category)

    ## return(list(call = call.dt,
    ##             contigs = final.tigs,
    ##             distal.bps = distal.bps,
    ##             all.bps = keep.tigs.breakends,
    ##             all.tigs = all.tigs))
## }

#' @name grab_distal_breakends
#' @title grab_distal_breakends
#'
#' @description
#'
#' Grab candidate distal breakends based on tumor-specificity
#'
#' If there are tumor-specific breakends:
#' If there are MAPQ = 60 breakends all within 1kbp: returns a single breakend (the first one)
#' If there are MAPQ = 60 breakends but not within 1 kbp: return MAPQ = 60 breakends
#' Otherwise, return all breakends
#'
#' @param bps (data.table) data.table representing breakpoints coercible to GRanges. bps needs to have columns mapq, tumor.count, and normal.count representing supporting barcodes for that contig.
#' @param mapq.thresh (default 60)
#' @param specificity.thresh (numeric) cutoff for tumor-specific contigs
#' @param pad (numeric) default 1 kbp
#'
#' @return a data.table with filtered breakends and additional metadata columns
grab_distal_breakends = function(bps, mapq.thresh = 60, specificity.thresh = 0.05, pad = 1e3) {

    if (!bps[, .N]) {
        return(bps)
    }

    ## check whether there are any tumor-specific breakends
    distal.bps = copy(bps)[(!proximal)]

    if (!distal.bps[, .N]) {
        return(distal.bps)
    }
    
    distal.bps[, clustered := NA]
    distal.bps[, specific := normal.count <= specificity.thresh * tumor.count & tumor.count > 1]

    if (!distal.bps[(specific), .N]) {
        return(distal.bps[(specific),])
    }

    distal.bps = distal.bps[(specific),]

    ## label high mapq breakends
    distal.bps[, high.mapq := mapq >= mapq.thresh]

    ## filter as much as possible
    if (distal.bps[(high.mapq), .N]) {
        distal.bps = distal.bps[(high.mapq)]
    }

    if (distal.bps[(proximal.first), .N]) {
        distal.bps = distal.bps[(proximal.first)]
    }
    
    ## prefer contigs with split read support
    if (distal.bps[(!mate.only), .N]) {
        distal.bps = distal.bps[(!mate.only)]
    }

    ## check if all distal breakends go to the same place
    bps.reduced.gr = GenomicRanges::reduce(dt2gr(distal.bps[, .(seqnames, start, end, strand)]) + pad) - pad
    if (length(bps.reduced.gr) == 1) {
        return(distal.bps[, clustered := TRUE])
    } else {
        return(distal.bps[, clustered := FALSE])
    }
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
#' @param bps (data.table) data.table representing breakpoints coercible to GRanges. bps needs to have columns mapq, tumor.count, and normal.count representing supporting barcodes for that contig.
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
#' - mapq (numeric)
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
        extra.cols = mcols( dt2gr(res)[, c()] %$% dt2gr(calns)[, c("c_type", "c_spec", "cigar", "mapq")] )
        res[, c_type := extra.cols$c_type]
        res[, c_spec := extra.cols$c_spec]
        res[, cigar := extra.cols$cigar]
        res[, mapq := extra.cols$mapq]
        res[, seed := gr.string(seed)]
    }
    return(res)
}

