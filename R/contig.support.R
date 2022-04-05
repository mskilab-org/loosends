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
#' @param min.supp.reads (numeric) min number of supporting reads to be considered as having read support (default 2)
#' @param max.normal.reads (numeric) maximum accepted number of supporting reads in normal (default 5)
#' @param min.read.ratio (numeric) minimum ratio of tumor to normal reads (default 5)
#' @param verbose (logical) default FALSE
check_contig_support_wrapper = function(calns,
                                        reads.dt,
                                        ref,
                                        name.field = "name",
                                        min.supp.reads = 2,
                                        max.normal.reads = 10,
                                        min.read.ratio = 10,
                                        verbose = FALSE)
{

    ## if empty datable, just return something with no rows
    if (!calns[, .N]) {
        return(data.table(name = character(),
                          tumor.count = numeric(),
                          normal.count = numeric(),
                          tumor.support = logical(),
                          tumor.specific = logical()))
    }

    ## get unique qnames
    if (!name.field %in% names(calns)) {
        stop("Invalid name field for contigs: ", name.field)
    }
    calns = copy(calns)[, name := get(name.field)]
    qns = unique(calns[, name])

    supporting.qname.counts = lapply(qns,
                                     function(qn) {
                                         if (verbose) { message("Checking qname: ", qn) }
                                         ## get seed frame 
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
                                             message("Number of tumor read pairs: ", tumor.count)
                                             message("Number of normal read pairs: ", normal.count)
                                         }

                                         tumor.support = tumor.count > min.supp.reads
                                         normal.support = normal.count > min.supp.reads
                                         tumor.specific = tumor.support
                                         if (normal.count > max.normal.reads) {
                                             tumor.specific = FALSE
                                         }
                                         if (normal.count > 0 & (tumor.count / normal.count < min.read.ratio)) {
                                             tumor.specific = FALSE
                                         }
                                         return(data.table(name = qn,
                                                           tumor.count = tumor.count,
                                                           normal.count = normal.count,
                                                           normal.support = normal.support,
                                                           tumor.support = tumor.support,
                                                           tumor.specific = tumor.specific))
                                     })

    supporting.qname.counts = setkey(rbindlist(supporting.qname.counts), "name")

    ## annotate calns with qname counts
    calns[, tumor.count := supporting.qname.counts[name, tumor.count]]
    calns[, normal.count := supporting.qname.counts[name, normal.count]]
    calns[, tumor.support := supporting.qname.counts[name, tumor.support]]
    calns[, normal.support := supporting.qname.counts[name, normal.support]]
    calns[, tumor.specific := supporting.qname.counts[name, tumor.specific]]
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
    
    rs = readsupport::contig.support(reads = window.reads.gr,
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
    
    rs = readsupport::contig.support(reads = window.reads.gr,
                                     contig = calns.gr,
                                     ref = ref,
                                     cg.contig = cg.contig,
                                     verbose = verbose)

    window.reads.dt[, supporting := qname %in% rs$qname]
    return(window.reads.dt[(supporting),])
}
