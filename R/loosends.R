#' @name junction.support
#' @title junction.support
#' @description
#'
#' Forked from skitools to allow manual input of reference contig
#'
#' Takes as input a GRanges of bam alignments (e.g. outputted from bamUtils::read.bam) and a GRanges of rearranged
#' reference aligned contigs (e.g. output of RSeqLib::BWA) and a set of Junction objects, and outputs reads supporting
#' these junctions by building a contig around each junction (from the reference) and then running contig.support (see
#' that functions docuemntation for criteria)
#'
#' @param loose.end.str (character)
#' @param reads GRanges in SAM / BAM format e.g. output of read.bam or BWA, with fields $qname, $cigar, $flag $seq all populated in standard fashion, and optionally $AS
#' @param junctions Junction object
#' @param bwa RSeqLib BWA object and path to fasta file corresponding to the reference
#' @param ref optional DNAStringSet corresponding to reference genome sequence
#' @param pad padding around the junction breakpoint around  which to analyze contig and reference sequences, this should be several standard deviations above the average insert size (2000)
#' @param realign flag whether to realign or just use existing alignments
#' @param bx logical flag whether data is linked reads, must then have BX flag, and the pad will be set to minimum 1e5
#' @param refseq.pad (1e3)
#' @param verbose logical flag (TRUE)
#' @param ... additional parameters to contig support
#' @return reads re-aligned to the reference through the contigs with additional metadata describing features of the alignment
#' @export
#' @author Marcin Imielinski
junction.support = function(loose.end.str = NA_character_, reads, junctions = NULL, bwa = NULL, ref = NULL, pad = 500, bx = FALSE, pad.ref = pad*20, refseq.pad = 1e3, realign = TRUE, walks = NULL, verbose = TRUE, ...)
{
    if (is.na(loose.end.str) | !inherits(loose.end.str, 'character')) {
        stop('must supply character string for loose end')
    }

    wholseed = parse.gr(loose.end.str) + pad

  if (!inherits(reads, 'GRanges') || is.null(reads$qname) || is.null(reads$cigar) || is.null(reads$seq) || is.null(reads$flag))
    stop('read input must be GRanges with fields $qname, $cigar, $seq, $flag and optionally $AS')

  sl = seqlengths(reads)
  if (bx)
    pad = max(pad, 1e5)

  if (!is.null(junctions))
    walks = jJ(junctions$grl)$gw(pad = pad)

  if (is.null(walks))
    stop('Either walks or junctions must be provided')

  if (bx)
  {
    if (is.null(reads$BX))
      stop('reads must have BX tag, may need to read.bam with tag option to extract it')

    if (!length(reads))
      return(reads)

    sc = score.walks(walks$grl, reads = reads, verbose = FALSE, raw = TRUE)$sc
    res = as.data.table(melt(as.matrix(sc)))[value>0, .(BX = Var1, walk = Var2)]
    reads = gr2dt(reads) %>% merge(res, by = 'BX') %>% dt2gr
    return(reads)
  }

  if (!realign)
  {
    if (is.null(junctions))
      junctions = walks$edges$junctions

    ## strand flip since 
    ## read orientation convention
    ## is opposite to junction convention
    reads = gr.flipstrand(reads) 
    reads$R1 = bamUtils::bamflag(reads$flag)[,'isFirstMateRead']>0
    r1 = reads %Q% (R1 == TRUE) %>% as.data.table
    r2 = reads %Q% (R1 == FALSE) %>% as.data.table
    ov = merge(r1, r2, by = 'qname')
    grl = grl.pivot(
      GRangesList(dt2gr(ov[, .(seqnames = seqnames.x, start = start.x, end =end.x, strand = strand.x)],
                        seqlengths = sl),
                  dt2gr(ov[, .(seqnames = seqnames.y, start = start.y, end = end.y, strand = strand.y)],
                        seqlengths = sl)))
    values(grl)$qname = ov$qname
    ## make junctions out of reads and cross with "real" junctions
    jn = merge.Junction(jJ(grl), junctions, cartesian = TRUE, pad = pad)
    if (!length(jn))
      return(reads[c()])
    out = merge(as.data.table(gr.flipstrand(reads)), unique(jn$dt[, .(qname, junction.id = subject.id)]), by = 'qname') %>% dt2gr(seqlengths = sl)
    return(out)
  }
  
  if (inherits(bwa, 'character') && file.exists(bwa))
  {
    if (verbose)
      message('Loading BWA index')
    bwa = BWA(bwa)
  }

  if (!inherits(ref, 'DNAStringSet'))
  {
    if (verbose)
      message('Loading genome reference as DNAStringSet')

    ref = rtracklayer::import(bwa@reference)
  }

  ## only use the fasta header before the first space as the seqnames of ref 
  names(ref) = strsplit(names(ref), '\\s+') %>% sapply('[', 1)

  if (length(setdiff(seqnames(walks$nodes$gr), seqlevels(ref))))
    stop('seqlevels mismatch between junctions / walks and reference, please check reference (e.g. chr issues)')

  if (length(setdiff(seqnames(walks$nodes$gr), seqlevels(bwa))))
    stop('seqlevels mismatch between junctions / walks and BWA reference, please check reference (e.g. chr issues)')

  if (verbose)
    message('Building and mapping derivative contigs')

  contig = bwa[ref[gr.fix(walks$grl, ref, drop = TRUE)]]

  if (verbose)
    message('Building reference contigs flanking loose end')
   ## contigref = ref[gr.fix(walks$edges$junctions$footprint + pad.ref, ref, drop = TRUE)]
    rs = getSeq(Hsapiens, trim(gr.fix(gr.chr(wholseed + refseq.pad), Hsapiens)))
    refseq = as.character(rs)


  if (verbose)
    message('Making gChain mapping contigs to reference')
  cg.contig = gChain::cgChain(contig)

  if (verbose)
    message('Running contig support')

    ## reads = contig.support(reads, contig, ref = contigref, cg.contig = cg.contig, ...)
    reads = contig.support(reads, contig, ref = refseq, cg.contig = cg.contig, ...)
  reads$junction.id = reads$contig.id
  return(reads)  
}

#' @name contig.support
#' @title contig.support
#' @description
#'
#' Takes as input a GRanges of bam alignments (e.g. outputted from bamUtils::read.bam) and a GRanges of rearranged
#' reference aligned contigs (e.g. output of RSeqLib::BWA).
#'
#' It identifies the subset of reads that support each of the contigs and "lifts" those reads
#' through the read --> contig and contig --> reference alignments, returning supporting reads in reference coordinates.
#'
#' The criteria for support include min.bases aligning to at least two chunks of the rearranged contig, and
#' requirement that min.aligned.frac fraction of bases in every supporting read is aligned to that contig.
#'
#' Additional requirements for support include not allowing split alignment of individual reads to the contigs
#' (note: this does not mean we don't detect split reads that support the structural variant, this is captured
#' by the contig -> reference alignment, we are just requiring the reads align (near) perfectly to the contig).
#' and requiring alla alignments from a read pair (oriented to R1 frame of the fragment) to align to the same
#' strand of the contig.
#'
#' Finally, reads are not included in support if they align better to the reference than their native alignment,
#' which is determined by comparing the $AS of their contig alignment with their original alignment score, stored
#' in the provided metadata $AS field.  If reference AS is not provided as metadata, it will is assumed to be zero. 
#'
#' $AS can be optionally recomputed against a DNAStringSet "ref" that represent the reference
#' sequence.  (Note that this "ref" does not have to be the full genome reference, it is just used to compute
#' the alignment scores, and in fact for this to work  efficiently, it's recommended that the provided
#' reference sequence is local to the regions of interest, e.g. a few kb flanking each SV breakpoint,
#' rather than the whole genome.)
#'
#' The outputted reads include additional metadata including number of bases aligning to each chunk of the aligned contig.
#' 
#' 
#' @param reads GRanges in SAM / BAM format e.g. output of read.bam or BWA, with fields $qname, $cigar, $flag $seq all populated in standard fashion, and optionally $AS
#' @param contig GRanges in SAM / BAM format wth fields $qname, $cigar and $seq all [populated
#' @param ref optional DNAStringSet representing a reference sequence to compute alignments against
#' @param chimeric logical flag whether to require reads to support junctions in chimericcontigs (ie discontiguous chunks on the reference), chimeric = FALSE
#' @param strict strict requires that the alignment score of the read to contig alignment needs to be better for at least one read (and also not worse for any of the reads) 
#' @param isize.diff (numeric) the insert size in the reference vs. the insert size in the contig between discordant read pairs
#' 
#' @return reads re-aligned to the reference through the contigs with additional metadata describing features of the alignment
#' @export
#' @author Marcin Imielinski

contig.support = function(reads, contig, ref = NULL, chimeric = TRUE, strict = TRUE, cg.contig = gChain::cgChain(contig), isize.diff = 1e3, min.bases = 20, min.aligned.frac = 0.95, new = TRUE, 
                          verbose = TRUE)
{
    if (length(reads)==0)
        stop('reads must be non empty GRanges with $qname, $cigar, $seq, and $flag fields')

    if (length(contig)==0)
        stop('contigs must be non empty GRanges with $qname, $cigar and $seq fields')

    if (verbose)
        message('Prepping reads for contig alignment')
    seq = unique(gr2dt(contig), by = c('qname'))[, structure(as.character(seq), names = as.character(qname))]
    bwa.contig = RSeqLib::BWA(seq = seq)
    chunks = gChain::links(cg.contig)$x
    strand(chunks) = '+'
    chunks = disjoin(chunks)
    if (is.null(reads$R1))
        reads$R1 = bamflag(reads$flag)[,'isFirstMateRead']>0
    reads$read.id = 1:length(reads)
    if (is.null(reads$AS))
    {
        warning('AS not provided in reads, may want to consider using tag = "AS" argument to read.bam or provide a ref sequence to provide additional specificity to the contig support')
        reads$AS = 0
    }
    nix = as.logical(strand(reads) == '-' )
    reads$seq[nix] = reverseComplement(DNAStringSet(reads$seq[nix])) ## flip read sequences to original strand
    reads[!reads$R1] = gr.flipstrand(reads[!reads$R1]) ## flip R2 read orientation to R1 strand
    reads$seq[!reads$R1] = reverseComplement(DNAStringSet(reads$seq[!reads$R1])) ## flip R2 read sequences to R1 strand
    reads = reads %Q% (!duplicated(paste(qname, R1)))

    if (!is.null(ref)) ## realign reads against reference DNAStringSet if provided to get alignment scores
    {
        if (verbose)
            message('Realigning reads against reference DNAStringSet')
        
        if (inherits(ref, 'character') | inherits(ref, 'DNAStringSet')) {
            if (verbose) {
                message("building BWA from sequence")
            }
            bwa.ref = RSeqLib::BWA(seq = ref)
        } else if (inherits(ref, 'BWA')) {
            if (verbose) {
                message("Using supplied BWA")
            }
            bwa.ref = ref
        } else {
            stop("invalid data type provided for ref: ", class(ref)[1])
        }

        tmp = bwa.ref[reads$seq] %>% gr2dt
        tmp$ix = as.numeric(as.character(tmp$qname))
        tmp$R1 = reads$R1[tmp$ix]
        tmp$qname = reads$qname[tmp$ix]
        tmp = unique(tmp, by = c('qname', 'R1'))
        setkeyv(tmp, c('qname', 'R1'))
        if (nrow(tmp))
        {
            ## what if not both read pairs are included?
            tmp[, isize := ifelse(any(seqnames != seqnames[1] | any(strand != strand[1])), NA_integer_, diff(range(start, end))), by = qname]
            ## isize should be NA if only one of the sides has a match
            unpaired.qnames = tmp[, .(count = .N), by = qname][count == 1, qname]
            tmp[qname %in% unpaired.qnames, isize := NA_integer_]
            ## the final isize should only account for the isize from the provided refseq
            reads$isize = pmin(tmp[.(reads$qname, reads$R1), isize], Inf, na.rm = TRUE)
            ## next line is the original version
            ## reads$isize = pmin(tmp[.(reads$qname, reads$R1), isize], reads$isize, Inf, na.rm = TRUE)
            ## reads[, c("isize", "isize.new", "isize.pmin", "sample")] %Q% (sample == "A0K8N")
            ## reads %Q% (qname %in% unpaired.qnames)
            reads$AS = tmp[.(reads$qname, reads$R1), AS]
            ## replace cigar
            reads$ref.aligned.cigar = tmp[.(reads$qname, reads$R1), cigar]
        }
    }

    if (verbose)
        message('Aligning reads against derivative contigs')
    

    ## aligning reads to contig
    rdt = as.data.table(reads)
    ## use new cigar?
    rdt[, ref.aligned := countCigar(ref.aligned.cigar)[, 'M']] 
    rdt[, ref.aligned := countCigar(cigar)[, 'M']] ## this is not the cigar from the realignment?
    rdt[, ref.aligned.frac := ref.aligned/qwidth[1], by = .(qname, R1)] 

    reads$ref.aligned.frac = rdt$ref.aligned.frac

    ## readsc is a data.table that denotes read locations in contig coordinates
    readsc = bwa.contig[reads$seq] %>% gr2dt
    readsc$cigar = as.character(readsc$cigar)
    readsc$ix = as.integer(as.character(readsc$qname))
    readsc$R1 = reads$R1[readsc$ix]
    readsc$read.id = reads$read.id[readsc$ix]

    ## these are splits on the contig, not reference --> shouldn't be any for good alignment
    readsc[, nsplit := .N, by = .(qname, R1)] 
    readsc[, aligned := countCigar(cigar)[, 'M']]

    ## these are splits on the contig, not reference --> shouldn't be any for good alignment
    readsc[, aligned.frac := aligned/qwidth[1], by = .(qname, R1)]
    readsc$AS.og = reads$AS[readsc$ix]
    readsc$isize = abs(reads$isize[readsc$ix])


    readsc$seqnames.og = seqnames(reads)[readsc$ix] %>% as.character
    readsc$strand.og = strand(reads)[readsc$ix] %>% as.character
    readsc$start.og = start(reads)[readsc$ix]
    readsc$end.og = end(reads)[readsc$ix]
    readsc$ref.isize = gr2dt(readsc)[, ref.isize := ifelse(
                                           all(seqnames.og == seqnames.og[1]) & all(strand.og == strand.og[1]),
                                           as.numeric(diff(range(c(start.og, end.og)))),                                   
                                           Inf), by = qname]$ref.isize %>% abs

    readsc$ref.aligned.frac = reads$ref.aligned.frac[readsc$ix]
    readsc$AS.og[is.na(readsc$AS.og)] = 0
    readsc$qname = reads$qname[readsc$ix]

    ## track sample (comment out later)
    tst = readsc[, .(aligned.frac, AS, AS.og, ref.aligned.frac, qname)][, sample := rdt$sample[match(qname, rdt$qname)]]
    tst[, .(aligned.frac, ref.aligned.frac, AS, AS.og, sample)]
    
    ## new scoring method based on cgChain of reads to contigs
    if (new)
    {
        ## cgChain representing read to contig alignments
        readsc$al.id = 1:nrow(readsc)

        if (verbose)
            message('Generating read to contig cgChain')
        alcg = gChain::cgChain(readsc)
        alchunks = cbind(as.data.table(values(alcg)),
                         as.data.table(gChain::links(alcg)$x),
                         as.data.table(gChain::links(alcg)$y)[, .(contig = seqnames,
                                                                  contig.start = start,
                                                                  contig.end = end,
                                                                  contig.strand = strand)])

        ## strands should be aligned to read / fragment + strand, but if not let's flip
        alchunks[strand == '-', ":="(strand = '+', contig.strand = c('+' = '-', '-' = '+')[contig.strand])]

        ## now for each al.id (ie bam record) let's pick the left most gChain / links record on the read / fragment 
        ## ie this is the lowest coordinate on the query
        ## (note that cgChain will split indels into separate ranges hence giving one to many mapping of al.id
        ## to records in links)
        setkeyv(alchunks, c('qname', 'start', 'end'))
        ## alchunks[, is.min := start == min(start), by = al.id]
        ## alchunks = alchunks[is.min == TRUE, ]

        ## so now we want to find alignments that are
        ## (1) concordant with respect to the contig
        ##  i.e. there is a monotonic increase (decrease) of contig.start if the contig.strand is + (-)
        ## (2) most of the read (aligned.frac) is represented
        ## (3) AS scores are better than original
        ## (4) isize better than original (where isize is the contig. span between the first and last alignment) .. related to (1)


        if (verbose)
            message('Scoring read to contig to alignments')
        alchunks[, contig.sign := ifelse(contig.strand == '+', 1, -1)]
        alchunks[, concordant.sign := all(contig.sign == contig.sign[1]), by = qname]

        ## check if we never go from R1 == FALSE to R1 == TRUE
        alchunks[, concordant.R1R2 := all(diff(!R1)>=0), by = qname]

        ## check to see that our contig.start always increasing or decreasing
        alchunks[, concordant.start := all(diff(contig.sign[1]*contig.start)>0), by = qname]

        alchunks[, contig.isize := diff(range(contig.start, contig.end)), by = qname]
        alchunks[, bases := sum(width), by = qname]

        ## nonzero width of chunks with better alignment scores relative to REF
        ## should this also be by read? (e.g. qname, R1)
        alchunks[, AS.better := sum(width[AS>AS.og]), by = .(qname, R1)]
        alchunks[, AS.worse := sum(width[AS<AS.og]), by = .(qname, R1)]
        alchunks[, AS.equal := sum(width[AS==AS.og]), by = .(qname, R1)]

        keepq = alchunks[concordant.sign & concordant.R1R2 & concordant.start &
                         bases > min.bases & aligned.frac > min.aligned.frac &
                         aligned.frac >= ref.aligned.frac &
                         nsplit  == 1 & ## no split reads
                         contig.isize - ref.isize < isize.diff & ## reimplementation of isize.diff
                         (AS.better>0 | contig.isize<ref.isize) & ## use isize to prevent removing fb's
                         AS.worse == 0, ] ## all bases non-inferior alignment to original

        ## keep read-specific information... 
        keepq = keepq[, .(qname, R1, contig, contig.isize, contig.strand, bases,
                          contig.sign, AS.better, AS.worse, AS.equal)] %>% unique(by = c('qname', 'R1'))
    }
    else ## old scoring method
    {
        ## if strict (default) remove any alignments that overlap others in the same qname
        if (strict)
        {
            readsc = dt2gr(readsc)
            readsc = readsc %Q% (rev(order(AS)))
            readsc = readsc[!duplicated(gr.match(readsc, readsc, by = 'qname')), ] %>% gr2dt
        }


        if (verbose)
            message('Computing overlap stats')

        ov = dt2gr(readsc) %*% chunks
        strand(ov) = readsc$strand[ov$query.id]
        ov$subject.id = paste0('chunk', ov$subject.id)
        ovagg = dcast.data.table(ov %>% gr2dt, qname ~ subject.id, value.var = 'width', fun.aggregate = sum)
        ovagg$nchunks = rowSums(ovagg[, -1]>min.bases)  ## good means we hit multiple chunks with sufficient bases
        rstats = gr2dt(ov)[, .(
                          contig.id = unique(seqnames)[1],
                          pos = sum(width[strand == '+']),
                          neg = sum(width[strand == '-']),
                          aligned.frac = min(aligned.frac),
                          num.contigs = length(unique(seqnames)), ### fixing later ... multiple contigs as input could distort results
                          paired = any(R1) & any(!R1), 
                          isize.contig = diff(range(c(start, end))),
                          isize.og = isize[1],
                          qsplit = any(nsplit>1), ## any sequences in this qname split on the contig ie a bad alignment on the contig
                          worse = any(AS.og>AS), ## any alignment in this qname worse than vs reference?
                          better = any(AS>AS.og) ## any alignment in this qname better than vs reference?
                      ), by = qname] %>% merge(ovagg, by = 'qname')

        ## apply filters ie nchunks>1 if chimeric, all alignments have to be of one sign
        ## if not paired then AS < AS.og else isize<isize.og
        keepq = rstats[nchunks>chimeric & (pos == 0 | neg  == 0) & aligned.frac > min.aligned.frac & !worse & (better | !strict | (paired & isize.contig < isize.og - isize.diff)) & !qsplit & num.contigs == 1, ]
        if (nrow(keepq)==0)
            return(reads[c()])

        keepq$aligned.frac = NULL
    }

    ## R1/R2-speicific information is lost here
    ##readsc = merge(readsc, keepq, by = 'qname') %>% dt2gr
    readsc = merge(readsc, keepq, by = c('qname', 'R1')) %>% dt2gr
    
    if (verbose)
        message('Lifting reads through contig back to reference')

    out = gChain::lift(cg.contig, readsc)

    if (length(out)) ## add reads metadata back to out
    {
        out[!out$R1] = gr.flipstrand(out[!out$R1])
        out$col = ifelse(out$R1, 'blue', 'gray')

        if (verbose)
            message('Adding metadata to reads')
        metacols = setdiff(names(values(reads)), names(values(out)))
        values(out) = cbind(values(out), values(reads)[match(out$read.id, reads$read.id), metacols])
    }

    if (verbose)
        message('Done')

    out
}

#' build.from.win
#'
#' assembles strand-specific reads and their mates into contigs and counts read support per contig within a single seed window
#' @param win GRanges seed window
#' @param ri data.table reads and their mates
#' @param tracks optional, character which tracks (sample + strand) to assemble, default all
#' @param align.thres (numeric) width of alignment / number of non-N characters
build.from.win = function(win, ri, tracks=NULL, verbose = FALSE, align.thres = 0.9){
    if(!("track" %in% colnames(ri))){
        ri[, track := ifelse(strand=="+", "for", "rev")]
    }
    if(is.null(tracks)) tracks = ri[, unique(track)]

    if (verbose) {
        message("Tracks used for assembly: ")
        message(paste(tracks, collapse = '\n'))
    }
    ri[track %in% tracks, {
        t = track[1]
        ri[, seed := dt2gr(ri) %N% win > 0 & track==t]
        rtmp = ri[qname %in% qname[seed]]
        rtmp[(seed)][is.dup(qname), seed := 1:.N == 1, by=qname] ## in case of a foldback where both reads are "seed", choose one as the anchor
        rtmp[, seed.frame := ifelse(seed, reading.frame, as.character(reverseComplement(DNAStringSet(reading.frame))))]
        if(nrow(rtmp) < 5){
            data.table(peak = gr.string(win[,c()]), seq = as.character(NA), good.assembly=FALSE, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads=nrow(rtmp))
        } else{
            srf = rtmp[, seed.frame]
            
            if (verbose) {
                message("Starting assembly with ", length(srf), " reads")
            }
            f = Fermi(srf, assemble=T)
            contigsf = contigs(f)
            unassembled.reads = length(srf)

            if (verbose) {
                message("Initial number of contigs: ", length(contigsf))
            }
            
            if (length(contigsf)) {
                reassemble = TRUE
                while(reassemble){

                    ## using current contigs, count the number of unassembled reads
                    ctigs = BWA(seq=contigsf)
                    aln = ctigs[srf]

                    ## qnames not included in alignment
                    s1 = which(!(seq_along(srf) %in% as.integer(aln$qname)))

                    ## qnames with "bad" alignments (e.g. alignment width < 0.9 * #non-N characters)
                    aln.dt = as.data.table(aln)
                    align.frac.dt = aln.dt[, .(best.aln = max(width / nchar(gsub("N", "", seq)))), by = qname]
                    s2 = align.frac.dt[best.aln < align.thres, as.integer(qname)]
                    ## s2 = as.data.table(aln)[, max(width / nchar(gsub("N", "", seq))), by=qname][V1 < 0.9, as.integer(qname)]
                    ## if the number of unassembled reads does not change from the previous iteration...
                    if(length(s1) + length(s2) == unassembled.reads | length(s1) + length(s2) < 7) {
                        reassemble = FALSE
                    } else {
                        unassembled.reads = length(s1) + length(s2)

                        if (verbose) {
                            message("Reassembling. Number of unassembled reads: ", unassembled.reads)
                        }
                        ## if(unassembled.reads > 7){
                        ## message("Starting assembly with ", unassembled.reads, " reads")
                        c2 = Fermi(srf[c(s1, s2)], assemble=T)
                        if(length(contigs(c2))){
                            contigsf = c(contigsf, contigs(c2))
                        } else {
                            reassemble = FALSE
                        }
                    }
                }
                goodassembly = rep(TRUE, length(contigsf))
                for(i in seq_along(contigsf)){
                    ctig = BWA(seq = contigsf[i])
                    ra = ctig[srf]
                    sc = table(factor(strand(ra), c("+", "-")))
                    if(!any(sc[1:2] == 0)) goodassembly[i] = FALSE
                    if(sc["+"] < sc["-"]) contigsf[i] = as.character(reverseComplement(DNAStringSet(contigsf[i])))
                }
                data.table(peak = gr.string(win[,c()]), seq = contigsf, good.assembly=goodassembly, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads=unassembled.reads)
            } else{
                data.table(peak = gr.string(win[,c()]), seq = as.character(NA), good.assembly=FALSE, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads = length(srf))
            }
        }
    }, by=track]
}

#' match.seq
#'
#' returns logical vector of length subject indicating whether matches were found with any query
#' @param query PDict of query sequences to match
#' @param subject DNAStringSet of sequences to parse for matches

match.seq = function(query, subject)
{
    if (is.null(names(subject)))
        names(subject) = 1:length(subject)

    if (any(duplicated(names(subject))))
    {
        warning('Names of subject sequences have duplicates, deduping')
        names(subject) = dedup(names(subject))
    }

    totals = lengths(vwhichPDict(query, subject)) > 0
    return(totals)
}

#' eighteenmer
#'
#' generate a PDict of telomeric 18mers
#' @param gorc 'g' returns 18mers of G-heavy strand only 'c' returns 18mers of C-heavy strand only 'both' returns 18mers of G-heavy strand and 18mers of C-heavy strand, default='both'
eighteenmer = function(gorc=c('g', 'c', 'both')){
    if(any(!(gorc %in% c('g', 'c', 'both')))) stop("Allowed values of gorc are 'g','c','both'")
    if(length(gorc)>1) gorc='both'
    pattern = system.file('extdata', 'telomeres.fa', package='loosends')
    query = tryCatch(readDNAStringSet(pattern), error = function(e) NULL)
    if (is.null(query))
    {
        query = DNAStringSet(readLines(pattern))
    }
    motifs = unname(as.character(query))
    if(gorc!='c'){
        g.motifs = motifs[grepl("GGG", motifs)]
        gquery = as.data.table(expand.grid(i=g.motifs,j=g.motifs,k=g.motifs,l=g.motifs,s=1:6))[, seq := substr(paste0(i, j, k,l), s, s+17)][, unique(seq)]
        gqueries = PDict(gquery)
        if(gorc=='g') return(gqueries)
    }
    if(gorc!='g'){
        c.motifs = motifs[grepl("CCC", motifs)]
        cquery = as.data.table(expand.grid(i=c.motifs,j=c.motifs,k=c.motifs,l=c.motifs,s=1:6))[, seq := substr(paste0(i, j, k,l), s, s+17)][, unique(seq)]
        cqueries = PDict(cquery)
        if(gorc=='c') return(cqueries)
    }
    return(PDict(c(gquery, cquery)))
}

#' munch
#' 
#' returns logical vector of length reads
#' indicating does read contain any match to query
#' searches for matches on same strand as input seq
#' @param reads GRanges or data.table containing $seq field with character sequence
#' @param query optional, PDict of motifs, default= 18mers of telomeric motifs (both strands)
#' @export
munch = function(reads, query=NULL){
    if(is.null(query)) query = eighteenmer()
    seq = reads$seq
    seq_frame_ref = DNAStringSet(seq)
    return(match.seq(query, seq_frame_ref))
}

#' filter.graph
#'
#' identifies quality loose ends from gGraph
#' @param gg gGraph of JaBbA genome graph (contains loose ends)
#' @param cov.rds character path to binned coverage data
#' @param purity optional, fractional purity of sample, default assumes 1
#' @param ploidy optional, ploidy of sample, default infers from gg
#' @param field optional, default = "ratio", column in binned coverage data
#' @param PTHRESH optional, threshold for GLM p-value for calling true positive loose ends, default=3.4e-7 provides consanguinity with large dataset bonferroni correction
#' @param verbose optional, default=FALSE
#' @return data.table containing a row for every loose end in gGraph gg with a logical column `true.pos` indicating whether each loose end has passed all filters (TRUE) or not (FALSE)
#' @export
filter.graph = function(gg, cov.rds, purity=NULL, ploidy=NULL, field="ratio", PTHRESH=3.4e-7, verbose=F){
    if(is.null(purity)) stop("must provide purity")
    ## gather loose ends from sample
    if(verbose) message("Identifying loose ends")
    ll = gr2dt(gr.start(gg$nodes[!is.na(cn) & loose.cn.left>0]$gr))[, ":="(lcn = loose.cn.left, strand = "+")]
    lr = gr2dt(gr.end(gg$nodes[!is.na(cn) & loose.cn.right>0]$gr))[, ":="(lcn = loose.cn.right, strand = "-")]
    l = rbind(ll, lr, use.names=T, fill=T)
    l[, leix := 1:.N]

    if(nrow(l)==0) {
        if(verbose) message("No loose ends")
        le.class = as.data.table(copy(l))[, true.pos := logical()]
    } else{
        ## call true positives
        l = dt2gr(l)
        if(verbose) message(sprintf("%s total loose ends in the graph", length(l)))
        le.class = filter.loose(gg, cov.rds, l, purity=purity, ploidy=ploidy, field=field, PTHRESH=PTHRESH, verbose=verbose)
    }
    return(le.class)
}

#' .waviness
#'
#' quantifies autocorrelation in coverage
.waviness = function(x, y, min.thresh = 5e3, max.thresh = 10e4, spar = 0.5, smooth = TRUE, filter = rep(FALSE, length(x)), trim = 10) {
    if(length(x)==0) return(NA)
    dat = data.table(x, y)[order(x), ]
    dat[, lag := x-min(x)]
    fdat = dat[!is.na(y), ][!is.infinite(y), ]
    ## autocorrelation
    if(nrow(fdat)==0) return(NA)
    fdat[, ac := as.numeric(acf(c(y, y), plot = FALSE, lag.max = length(y))$acf[-1])]
    if (smooth) { ## smoothing the autocorrelation gets rid of some more noise
        fdat = fdat[!is.na(lag) & !is.na(ac),]
        if (nrow(fdat[!is.na(lag) & !is.na(ac),])<4)
            return(NA)
        fdat$ac = predict(smooth.spline(fdat$lag, fdat$ac, spar = spar), fdat$lag)$y
    }
    return(fdat[lag>min.thresh & lag<max.thresh, sum(ac^2)])
}

#' .mod
#'
#' fit data table with glm and return residuals
.mod = function(dt){
    mod = dt[, glm(counts ~ tumor + fused + ix, family='gaussian')]
    res = dt$counts - predict(mod, dt, type='response')
    return(res)
}

#' .mod2
#'
#' fit data table with glm and return residuals
.mod2 = function(dt){
    mod = dt[, glm(counts ~ tumor + ix, family='gaussian')]
    res = dt$counts - predict(mod, dt, type='response')
    return(res)
}

#' filter.loose
#'
#' analyze coverage surrounding given loose ends to evaluate quality
#' @param gg gGraph of JaBbA model
#' @param cov.rds character path to binned coverage data
#' @param l data.table of loose ends to evaluate
#' @param jabba_rds (character) load gg and l directly from jabba output
#' @param purity optional, fractional purity of sample, default infers from gg
#' @param ploidy optional, ploidy of sample, default infers from gg
#' @param field optional, column name in cov.rds, default="ratio"
#' @param PTHRESH optional, threshold for GLM p-value for calling true positive loose ends, default=3.4e-7 provides consanguinity with large dataset bonferroni correction
#' @param waviness (logical) apply waviness filter? default TRUE
#' @param cov.delta (logical) apply filters on coverage ratio delta? default TRUE
#' @param tumor.delta (logical) apply filters on coverage ratio delta? default TRUE
#' @param normal.delta (logical) apply filters on coverage ratio delta? default TRUE
#' @param fdr (numeric) false discovery rate, default 0.05
#' @param epgap (numeric) epgap threshold, default 1e-3
#' @param ascat (character) path to ASCAT segmentation, default /dev/null
#' @param verbose optional, default=FALSE
#' @return data.table containing a row for every input loose end and logical column `true.pos` indicating whether each loose end has passed all filters (TRUE) or not (FALSE)
#' @export
filter.loose = function(gg,
                        cov.rds,
                        l,
                        jabba_rds = NULL,
                        purity=NULL,
                        ploidy=NULL,
                        field="ratio",
                        PTHRESH=3.4e-7,
                        waviness = TRUE,
                        cov.delta = TRUE,
                        tumor.delta = TRUE,
                        normal.delta = TRUE,
                        p.filter = TRUE,
                        fdr = 0.05,
                        epgap = 0.001,
                        ascat = "/dev/null",
                        verbose=F){
    ## load coverage and beta (coverage CN fit)
    if(verbose) message("Loading coverage bins")

    stopifnot(file.exists(cov.rds))

    cov = readRDS(cov.rds)
    cov = gr.sub(cov, "chr", "")
    
    if(!(field %in% colnames(values(cov)))) stop("must provide field in cov.rds")
    if(field != "ratio") cov$ratio = values(cov)[, field]
    if(!("tum.counts" %in% colnames(values(cov)))){
        yf = ifelse("reads.corrected" %in% colnames(values(cov)), "reads.corrected", field)
        cov$tum.counts = values(cov)[, yf]
    }
    if(!("norm.counts" %in% colnames(values(cov)))){
        cov$norm.counts = 1 ## dummy to make it flat
    }
    if(is.null(purity)) purity = gg$meta$purity ##stop("must provide purity") ## purity = 1
    if(is.null(ploidy)) ploidy = weighted.mean(gg$nodes$gr$cn, gg$nodes$gr %>% width, na.rm=T)
    ratios = cov$ratio
    beta = mean(ratios[is.finite(ratios)], na.rm=T) * purity/(2*(1-purity) + purity * ploidy)
    segs = gg$nodes$gr
    segs = gr.sub(segs, "chr", "")
    l = gr.sub(l, "chr", "")

    ## add leix to loose ends if not already there
    if (is.null(l$leix)) {
        l$leix = 1:length(l)
    }


    ## identify nodes flanking each loose end, extending up to 100kb away
    o = gr.findoverlaps(segs, l+1)
    segs = segs[o$query.id]; segs$leix = l[o$subject.id]$leix
    sides = gr.findoverlaps(segs, l+1e5, by="leix")
    values(sides) = cbind(values(sides), values(l[sides$subject.id]))
    sides$fused = !is.na(gr.match(sides, l, by="leix"))
    sides$wid = width(sides)
    sides$node.id = segs[sides$query.id]$node.id
    sides$cn = segs[sides$query.id]$cn

    ## gather coverage bins corresponding to fused & unfused sides of loose ends
    if(verbose) message("Overlapping coverage with loose end fused and unfused sides")
    rel = gr.findoverlaps(cov, sides)
    values(rel) = cbind(values(cov[rel$query.id]), values(sides[rel$subject.id]))
    qq = 0.05
    rel = gr2dt(rel)[, ":="(
                   in.quant.r = is.finite(ratio) & ratio >= quantile(ratio, qq, na.rm=T) & ratio <= quantile(ratio, 1-qq, na.rm=T),
                   good.cov=sum(is.na(tum.counts))/.N < 0.1 & sum(is.na(norm.counts))/.N < 0.1 & sum(is.na(ratio))/.N < 0.1 & wid > 5e4
               ), by=.(subject.id, fused)]

    rel[, lxxx := leix]
    variances = rel[(in.quant.r), var(ratio), keyby=.(fused, lxxx)]
    variances[, side := ifelse(fused, "f_std", "u_std")]
    variances[, std := sqrt(V1)]
    vars = dcast.data.table(variances, lxxx ~ side, value.var="std")
    rel[is.na(in.quant.r), in.quant.r := FALSE]
    rel[, tum.median := median(tum.counts[in.quant.r]), by=.(lxxx)]
    rel[, norm.median := median(norm.counts[in.quant.r]), by=.(lxxx)]
    rel[, tum.res := tum.counts - tum.median]
    rel[, norm.res := norm.counts - norm.median]
    tum.ks = rel[(in.quant.r), tryCatch(dflm(ks.test(tum.res[fused], tum.res[!fused])), error = function(e) dflm(ks.test(tum.res, tum.res))), by=lxxx][, p, by=lxxx]
    norm.ks = rel[(in.quant.r), tryCatch(dflm(ks.test(norm.res[fused], norm.res[!fused])), error = function(e) dflm(ks.test(norm.res, norm.res))), by=lxxx][, p, by=lxxx]
    pt1 = merge(vars, merge(tum.ks, norm.ks, by="lxxx", suffixes=c("_tum", "_norm"),all=T), by="lxxx", all=T)
    pt1$lxxx = as.character(pt1$lxxx); setkey(pt1, lxxx)
    pt1[, n_fdr := p.adjust(p_norm, "bonferroni")]
    pt1[, t_fdr := p.adjust(p_tum, "bonferroni")]

    rel[, ":="(
        tumor.mean.fused = mean(tum.counts[fused], na.rm=T),
        tumor.mean.unfused = mean(tum.counts[!fused], na.rm=T),
        normal.mean.fused = mean(norm.counts[fused], na.rm=T),
        normal.mean.unfused = mean(norm.counts[!fused], na.rm=T)
    ), by=leix]

    ## evaluate waviness across bins per loose end
    rel[, good.cov := all(good.cov), by=subject.id]
    if(verbose) message("Calculating waviness around loose end")
    rel[, waviness := max(.waviness(start[fused], ratio[fused]), .waviness(start[!fused], ratio[!fused]), na.rm=T), by=subject.id]

    ## prep glm input matrix
    if(verbose) message("Prepping GLM input matrix")
    glm.in = melt.data.table(rel[(in.quant.r),], id.vars=c("leix", "fused"), measure.vars=c("tum.counts", "norm.counts"), value.name="counts")[, tumor := variable=="tum.counts"]
    glm.in[, ix := 1:.N, by=leix]
    rel2 = copy(glm.in)
    setnames(glm.in, "leix", "leix2")

    ## calculate residuals from glm 
    rel2[, residual := .mod(glm.in[leix2==leix[1],]), by=leix]

    ## evaluate KS-test on residuals and calculate effect size
    ## effect will be from KS test on residuals
    ## estimate is replaced with difference of median coverage
    if(verbose) message("Running KS-test on fused vs unfused sides of loose ends")
    res = rel2[(tumor), tryCatch(dflm(ks.test(residual[fused], residual[!fused])), error=function(e) dflm(ks.test(residual, residual))), by=leix]
    est = rel2[, median(counts), keyby=.(fused, tumor, leix)][, V1[tumor] / V1[!tumor], keyby=.(leix, fused)][, V1[fused]-V1[!fused], keyby=leix]
    res$leix = as.character(res$leix); setkey(res, leix)
    res[as.character(est$leix), estimate := est$V1]

    ## combine relevant fields from each test
    test = rel2[(tumor), mean(counts), keyby=.(fused, leix)][, V1[fused]-V1[!fused], keyby=leix]; setnames(test, "V1", "testimate")
    test$leix = as.character(test$leix); setkey(test, leix)
    nest = rel2[!(tumor), mean(counts), keyby=.(fused, leix)][, V1[fused]-V1[!fused], keyby=leix]; setnames(nest, "V1", "nestimate")
    nest$leix = as.character(nest$leix); setkey(nest, leix)
    cnl = rev(rev(colnames(rel))[1:6])
    cns = c("name", "method", "estimate", "effect", "p")
    le.class = cbind(gr2dt(l[rel[!duplicated(subject.id), subject.id]]), rel[!duplicated(subject.id), cnl, with=F], res[rel[!duplicated(subject.id), as.character(leix)], cns, with=F], nest[rel[!duplicated(subject.id), as.character(leix)], "nestimate"], test[rel[!duplicated(subject.id), as.character(leix)], "testimate"])[, effect.thresh := beta]

    le.class[, f.std := pt1[.(as.character(leix)), f_std]]
    le.class[, u.std := pt1[.(as.character(leix)), u_std]]
    le.class[, ":="(n_fdr = pt1[.(as.character(leix)), n_fdr], t_fdr = pt1[.(as.character(leix)), t_fdr])]
    le.class[, bon := p.adjust(p, "bonferroni")]

    ## add JaBbA CN of fused and unfused sides
    fused.sides = as.data.table(sides)[fused == TRUE,]
    unfused.sides = as.data.table(sides)[fused == FALSE,]
    le.class[, fused.cn := fused.sides$cn[match(leix, fused.sides$leix)]]
    le.class[, unfused.cn := unfused.sides$cn[match(leix, unfused.sides$leix)]]

    ## correct p values
    if(verbose) message("Identifying true positives")

    ## check which loose ends overlap with ASCAT
    if ("passed" %in% colnames(le.class)) {
        le.class$passed = NULL
    }

    if ("true.pos" %in% colnames(le.class)) {
        le.class$true.pos = NULL
    }

    le.class[, passed := TRUE]

    ascat.le = numeric()
    
    if (file.exists(ascat) & file.info(ascat)$size) {
        if (verbose) {
            message("Checking overlap with ASCAT input")
        }
        ascat.gr = rtracklayer::import(ascat)
        ascat.bp = c(gr.start(ascat.gr), gr.end(ascat.gr))
        ascat.ov = gr.findoverlaps(gr.stripstrand(l)+1e3, ascat.bp, ignore.strand = TRUE)
        if (length(ascat.ov)) {
            ascat.le = unique(l[ascat.ov$query.id])$leix
        }
        le.class[!(leix %in% ascat.le), passed := FALSE]
    } else {
        if (verbose) {
            message("No ASCAT segmentation provided, skipping.")
        }
    }

    le.class[tumor.mean.fused < tumor.mean.unfused, passed := FALSE]
    le.class[fused.cn <= unfused.cn, passed := FALSE]

    if (waviness) {
        le.class[waviness >= 2, passed := FALSE]
    }

    if (p.filter) {
        le.class[is.na(p), passed := FALSE]
        le.class[p > PTHRESH, passed := FALSE]
    }
    if (cov.delta) {
        le.class[estimate < 0.6 * effect.thresh, passed := FALSE]
        le.class[estimate < f.std, passed := FALSE]
        le.class[estimate < u.std, passed := FALSE]
    }

    if (tumor.delta) {
        le.class[testimate < 0.6 * effect.thresh, passed := FALSE]
        le.class[t_fdr > fdr, passed := FALSE]
    }

    if (normal.delta) {
        le.class[abs(nestimate) > 0.6 * effect.thresh, passed := FALSE]
        le.class[n_fdr < fdr, passed := FALSE]
    }

    le.class[, true.pos := passed]


    ## if(!("epgap" %in% colnames(le.class))){
    ##     le.class[, passed := !is.na(p) & p < PTHRESH & estimate > (0.6*effect.thresh) & testimate > (0.6*effect.thresh) & waviness < 2 & abs(nestimate) < (0.6*effect.thresh)]
    ## }else le.class[, passed := !is.na(p) & p < PTHRESH & estimate > (0.6*effect.thresh) & testimate > (0.6*effect.thresh) & waviness < 2 & abs(nestimate) < (0.6*effect.thresh) & epgap < 1e-3]
    ## le.class[, true.pos := passed & estimate > f.std & estimate > u.std & n_fdr > 0.05 & t_fdr < 0.01]
    ## le.class$passed = NULL
    
    gc()
    return(le.class)
}

#' return.dt
#'
#' constructs output data.table from logical arguments
return.dt = function(reference, complex, missedj, novel, mystery, insertion, refmap=NULL, novmap=NULL, rrep=NULL, irep=NULL, nrep=NULL, junction=""){
    if(length(rrep)) if(length(rrep %Q% (c_spec == "unassembled"))) rrep[rrep$c_spec=="unassembled"]$c_spec = "low MAPQ"
    if(length(nrep)) if(length(nrep %Q% (c_spec == "unassembled"))) nrep[nrep$c_spec=="unassembled"]$c_spec = "low MAPQ"
    if(length(irep)) if(length(irep %Q% (c_spec == "unassembled"))) irep[irep$c_spec=="unassembled"]$c_spec = "low MAPQ"
    if(length(rrep)) if(length(rrep %Q% (c_spec != "low MAPQ"))) rrep = rrep %Q% (c_spec != "low MAPQ")
    if(length(nrep)) if(length(nrep %Q% (c_spec != "low MAPQ"))) nrep = nrep %Q% (c_spec != "low MAPQ")
    if(length(irep)) if(length(irep %Q% (c_spec != "low MAPQ"))) irep = irep %Q% (c_spec != "low MAPQ")

    frep = c(irep, nrep)
    rreps = ifelse(length(rrep), (rrep %Q% (rev(order(width))))[1]$c_spec, "")
    ireps = ifelse(length(irep), (irep %Q% (rev(order(width))))[1]$c_spec, "")
    nreps = ifelse(length(nrep), (nrep %Q% (rev(order(width))))[1]$c_spec, "")
    freps = ifelse(length(frep), (frep %Q% (rev(order(width))))[1]$c_spec, "")
    if(mystery){
        nreps = freps = "mystery"
        novmap=FALSE
    }
    call.string = NULL
    if(missedj) call.string = "missed junction"
    if(complex) call.string = paste(c(call.string, "complex"), collapse="; ")
    if(reference){
        radd = paste(ifelse(refmap, "mappable", "unmappable"), "seed repeat")
        if(rreps == "") rreps = "low MAPQ"
        radd = paste0(radd, ":", rreps)
        call.string = paste(c(call.string, radd), collapse="; ")
    }
    if(novel){
        nadd = paste(ifelse(novmap, "mappable", "unmappable"), "mate repeat")
        if(nreps == "") nreps = "low MAPQ"
        nadd = paste0(nadd, ":", nreps)
        call.string = paste(c(call.string, nadd), collapse="; ")
    }
    if(insertion){
        iadd = "insertion"
        if(ireps != "") iadd = paste0(iadd, ":", ireps)
        call.string = paste(c(call.string, iadd), collapse="; ")
    }
    if(mystery) call.string = paste(c(call.string, "mystery"), collapse="; ")
    if(is.null(refmap)) refmap = logical()
    if(is.null(novmap)) novmap = logical()
    ##    single.call = ifelse(complex, "complex rearrangement", ifelse(missedj, "missed junction", ifelse(reference & !(refmap), "seed repeat", ifelse(novel & !(novmap), "mate repeat", "mystery"))))
    single.call = ifelse(complex, "complex rearrangement", ifelse(missedj, "missed junction", ifelse(mystery | grepl("mystery", freps), "mystery", paste("type", as.integer(!refmap) + as.integer(!novmap), "loose end"))))
    if(single.call == "mystery") mystery = TRUE
    data.table(seedrep=reference,
               complex=complex,
               missedj=missedj,
               materep=novel,
               mystery=mystery,
               insertion=insertion,
               seed.repeats = rreps,
               ins.repeats = ireps,
               nov.repeats = nreps,
               mate.repeats = freps,
               seed.mappable = refmap,
               mate.mappable = novmap,
               junction = junction,
               call = call.string,
               category = single.call)
}

#' update.call
#'
#' modifies call string if logical column values have changed
#' (used when supplementing contig calls with discordant read pairs)
update.call = function(dt){
    dt = copy(dt)
    dt[, i := 1:.N]
    call.string = dt[, {
        call.string = character()
        if(missedj) call.string = "missed junction"
        if(complex) call.string = paste(c(call.string, "complex"), collapse="; ")
        if(seedrep){
            radd = paste(ifelse(seed.mappable, "mappable", "unmappable"), "seed repeat")
            radd = paste0(radd, ":", seed.repeats)
            call.string = paste(c(call.string, radd), collapse="; ")
        }
        if(materep){
            nadd = paste(ifelse(mate.mappable, "mappable", "unmappable"), "mate repeat")
            nadd = paste0(nadd, ":", nov.repeats)
            call.string = paste(c(call.string, nadd), collapse="; ")
        }
        if(insertion){
            iadd = "insertion"
            if(ins.repeats != "") iadd = paste0(iadd, ":", ins.repeats)
            call.string = paste(c(call.string, iadd), collapse="; ")
        }
        if(mystery) call.string = paste(c(call.string, "mystery"), collapse="; ")
        call.string
    }, by=i]
    return(call.string$V1)
}

#' merge.telomeres
#'
#' construct a GRanges in contig coordinates for telomere motif matches
merge.telomeres = function(treps, x, ch, out.dt, c_spec = "telomere"){
    values(treps) = out.dt[.(as.character(seqnames(treps))), colnames(ch@values), with=F]
    treps$AS = width(treps)
    treps$y = gr.string(dt2gr(as.data.table(treps)[, seqnames := "tel"]))
    treps$cigar = ifelse(start(treps) > 1, paste0(start(treps)-1, "S", width(treps), "M"), paste0(width(treps), "M"))
    treps$c_type = "rep"
    treps$mapq = 0
    treps$c_spec = c_spec
    treps$seed = FALSE
    x = dt2gr(rbind(as.data.table(x), as.data.table(treps), fill=TRUE, use.names=TRUE))
    return(x)
}

#' find.telomeres
#'
#' search for telomeric motif matches in contigs
find.telomeres = function(query, out.dt){
    gd = munch(out.dt, query)
    if(!(any(gd))) return(NULL)
    qns = out.dt[gd, qname]
    do.call('c', lapply(qns, function(qn) reduce(GRanges(qn, IRanges(do.call('c', lapply(vwhichPDict(query, DNAStringSet(out.dt[.(qn), seq]))[[1]], function(i) as.integer(gregexpr(as.character(query@dict0)[i], out.dt[.(qn), seq])[[1]]))), width=18), strand="+"))))
}

#' posterity.caller
#'
#' saving in case I mess up the new one:
#' parses contig alignments to assign loose end to category and describe repeat types
#' @param li data.table loose end to evaluate
#' @param calns optional, data.table of contig alignments, default=NULL (will not parse alignments)
#' @param insert optional, integer pad representing insert length in bp to identify seed alignments, default=750
#' @param pad optional, window around loose end to allow contig seed windows, default=1e3
#' @param uannot optional, GRanges of unmappable annotations, default bin/101.unmappable.annotations.rds
#' @param ref_dir optional, path to directory of unzipped reference tarball, default assumes 'package/extdata/hg19_looseends'
#' @param ref_obj optional, list of BWA objects built from ref_dir fastas, names must match expected "human" "rep" "polyA" "microbe" (only "human" is used), default=NULL
posterity.caller = function(li, calns = NULL, insert = 750, pad=NULL, uannot=NULL, ref_dir=system.file('extdata', 'hg19_looseends', package='loosends'), ref_obj=NULL) {
    reference = FALSE
    complex = FALSE
    missedj = FALSE
    novel = FALSE
    mystery = FALSE
    insertion = FALSE
    refmap = FALSE
    novmap = FALSE
    if(is.null(pad)) pad = 1e3
    if(is.null(uannot)){
        uannot = readRDS(system.file('extdata', '101.unmappable.annotations.rds', package='loosends'))
    }
    if(!is.null(ref_obj)) { human = ref_obj$human
    } else{
        if(!file.exists(ref_dir)) stop("Provide correct ref_dir containing reference .fa files")
        if(!file.exists(paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37_decoy.fasta")
        human = BWA(fasta=paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))
    }
    
    if(is.null(calns)) {
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            rrep = uannot %&% dt2gr(li)
            if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
        } else {
            refmap = T
            rrep = uannot[0]
        }
        mystery = mystery | !any(missedj, complex, reference, novel)
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, refmap=refmap, novmap=novmap, rrep=rrep))
    }
    if(!nrow(calns)){
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            rrep = uannot %&% dt2gr(li)
            if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
        } else {
            refmap = TRUE
            rrep = uannot[0]
        }
        mystery = mystery | !any(missedj, complex, reference, novel)                
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, refmap=refmap, novmap=novmap, rrep=rrep))
    }
    if(inherits(calns$mapq, "character")) calns$mapq = as.integer(calns$mapq)
    calns[c_spec == "polyA", c_type := "rep"]
    calns[c_spec == "ribosomal", c_type := "ribosomal"]
    seed = GRanges(calns$peak)
    strand(seed) = ifelse(grepl("for", calns$track), "+", "-")
    seed$qname = calns$qname
    ## narrow down to 1kb window on correct strand
    ## this filters out some that were only caught on the reverse ....
    calns = calns[!is.na(gr.match(seed, gr.flipstrand(dt2gr(li))+pad, ignore.strand=F))]
    if(!nrow(calns)){
        rrep = NULL
        if(dt2gr(li) %N% uannot) {
            reference = TRUE
            rrep = uannot %&% dt2gr(li)
            if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
        }
        mystery = mystery | !any(missedj, complex, reference, novel)
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, refmap=refmap, novmap=novmap))
    }

    ## filter out contigs with too few good input reads,
    ## unless there wouldn't be any left!
    ## keep going regardless!
    if(calns[map60.cov > 9, .N == 0]){
        reference = TRUE
        rrep = GRanges(li$seqnames, IRanges(li$start, li$end), c_spec="low MAPQ")
    }

    ## have added empty contigs for loose ends with none, to check the mapq of input reads
    calns = calns[!is.na(seqnames)]
    if(!nrow(calns)){
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            rrep = uannot %&% dt2gr(li)
            if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
        }
        mystery = mystery | !any(missedj, complex, reference, novel)
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, refmap=refmap, novmap=novmap))
    }
    
    ## checking this isn't just the seed region fully mapping & not going anywhere
    ## also making sure it doesn't leave but come back
    seed = GRanges(calns$peak)
    strand(seed) = ifelse(grepl("for", calns$track), "+", "-")
    seed$qname = calns$qname
    calns$seed = !is.na(gr.match(dt2gr(calns), seed+insert, by="qname", ignore.strand=F))
    if(any(calns$seed) & all(calns[(seed), mapq == 60])) refmap = TRUE
    ch = cgChain(calns)
    x = ch@.galx
    values(x) = ch@values
    x$y = gr.string(dt2gr(calns)[floor(as.numeric(rownames(ch@values)))])
    ## ignoring indel cigar alignments -- w/r/t when width is used to prioritize repeats
    x$source.id = floor(as.numeric(rownames(ch@values)))
    x = dt2gr(as.data.table(x)[, ":="(start = min(start), end = max(end)), by=source.id][!duplicated(source.id)][,!"source.id"])
    
    ## now that we're in contig coordinates, parse for telomeres
    out.dt = calns[!duplicated(qname)]
    setkey(out.dt, qname)
    tgreps = find.telomeres(eighteenmer('g'), out.dt)
    tcreps = find.telomeres(eighteenmer('c'), out.dt)
    if(!is.null(tgreps)) x = merge.telomeres(tgreps, x, ch, out.dt, c_spec="G telomere")
    if(!is.null(tcreps)) x = merge.telomeres(tcreps, x, ch, out.dt, c_spec="C telomere")
    ## viral calls overlapping telomere matches are ignored -- viral reference seqs include (TTAGGG)n
    if(any(x$c_type=="viral") & (!is.null(tgreps) | !is.null(tcreps))){
        virs = x %Q% (c_type == "viral")
        kvir = virs %O% (x %Q% (grepl("telomere", c_spec)))
        x = c(x %Q% (c_type != "viral"), virs[kvir < 0.9])
    }
    ## ribosomal calls overlapping polyA calls are ignored -- ribosomal ref includes (A)n
    if(any(x$c_spec=="ribosomal") & any((x %Q% (c_spec=="ribosomal")) %N% (x %Q% (c_spec=="polyA")))){
        ribs = x %Q% (c_spec == "ribosomal")
        krib = ribs %O% (x %Q% (c_spec == "polyA"))
        x = c(x %Q% (c_spec != "ribosomal"), ribs[krib < 0.9])
    }
    ## prioritze telomere motif matches over telomeric coordinate alignments
    if(any(x$c_spec=="Telomere / sub") & (!is.null(tgreps) | !is.null(tcreps))){
        tela = x %Q% (c_spec == "Telomere / sub")
        utel = is.na(gr.match(tela, x %Q% (grepl("tel", y)), by="seed"))
        x = c(x %Q% (c_spec != "Telomere / sub"), tela[utel])
    }

    
    ## make sure the contig goes and stays away:
    ## last 20 bases of the contig cannot align to seed window
    ends = gr.end(si2gr(x), width=20)
    ends$seed = ends %N% (x %Q% (seed))
    ends$human60 = ends %N% (x %Q% (c_spec=="human") %Q% (mapq==60) %Q% (!(seed)))
    ## outp: qnames of contigs with a seed alignment in the last 20 bases
    outp = as.data.table(ends)[seed > 0, as.character(seqnames)]

    if(reference){
        rrep = copy(x[1])
        rrep$c_spec = "low MAPQ"
    } else{
        rrep = x[0]
    }

    ## add to reference repeats: contig seed regions that also align elsewhere, at least 20bp
    ## (should be making sure the overlap specifically is at least 20bp)
    rrep = c(rrep, x %&% (x %Q% (seed)) %Q% (!(seed)) %Q% (width > 19) %Q% (seqnames %in% outp))
    if(any(rrep$c_spec == "human" & rrep$mapq < 60))
        rrep[rrep$c_spec == "human" & rrep$mapq < 60]$c_spec = "low MAPQ"
    if(length(rrep)) {
        rrep = rrep %Q% (c_spec!="human") %Q% (c_type!="human" | all(c_type=="human")) %Q% (!duplicated(c_spec))
        reference = TRUE
    }

    ## filter out contigs with seed alignment in the last 20 bases
    x2 = x %Q% (!(seqnames %in% outp)) ## save x so you can check for reference repeats
    if(!length(x2)){
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            sup = uannot %&% dt2gr(li)
            if(is.na(sup$c_spec)) sup$c_spec = sup$repClass
            rrep = grbind(rrep, sup)
        } else refmap = refmap | all((x %Q% (seed))$mapq == 60)
        mystery = mystery | !any(missedj, complex, reference, novel)
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, refmap=refmap, novmap=novmap))
    } else x = x2
    

    ## remember you won't always have caught the junction -- might just be the other side
    sections = disjoin(x)
    sections$score = sections %N% x
    sections$seeds = sections %N% (x %Q% (seed))
    sections$hg19 = sections %N% (x %Q% (c_type=="human"))
    sections$human = sections %N% (x %Q% (c_spec=="human"))
    sections$human60 = sections %N% (x %Q% (c_spec=="human") %Q% (mapq==60))

    ## SHOULD DO INSTEAD: trim the lowMQ bit to be only outside the repeat
    if(any(x$c_type=="rep") & any((x %Q% (c_type=="rep")) %N% (x %Q% (c_spec=="human") %Q% (mapq < 60)))){
        mq = x %Q% (c_spec == "human") %Q% (mapq < 60)
        kmq = mq %O% (x %Q% (c_type == "rep"))
        x = c(x %Q% (c_spec != "human" | mapq==60), mq[kmq < 0.9])
    }

    x = x %*% sections
    x = as.data.table(x)
    if(x[, any(seed)] & x[(seed), all(mapq<60)]){
        reference = TRUE
        rrep = grbind(rrep, dt2gr(copy(x[(seed)][1])[, c_spec := "low MAPQ"]))
    }
    x = dt2gr(x)
    sections = sections %$% (x %Q% (c_spec != "human"))[, 'c_spec']
    sections = dt2gr(as.data.table(sections))

    nrep = x[0]
    irep = x[0]

    ## can only catch a reference repeat this way if the seed is part of the contig
    rr = unique(as.character(seqnames(x %Q% (seeds>0) %Q% (score > 1) %Q% (width > 19))))
    if(length(rr)){
        reference = TRUE
        rrep = x %Q% (seqnames %in% rr) %Q% (seeds>0) %Q% (!(seed) | seeds > 1) %Q% (width > 19)
        if(any(rrep$c_spec == "human" & rrep$mapq < 60))
            rrep[rrep$c_spec == "human" & rrep$mapq < 60]$c_spec = "low MAPQ"
        if(length(rrep)){
            wides = as.data.table(gr.reduce(rrep, by="query.id"))[, query.id[width > 19]]
            rrep = rrep %Q% (c_spec!="human") %Q% (c_type!="human" | hg19==score) %Q% (query.id %in% wides) %Q% (!duplicated(c_spec))
        }
    }    

    ## make sure if there are multiple they go to the same place..
    ## seed has good MAPQ -- if not included in contig, align using getSeq
    ## last 20 bases overlap good MAPQ
    ## no MAPQ60 gaps 100bp+
    ## first check gaps:
    mj = setdiff(levels(seqnames(x)), unique(as.character(seqnames(gaps(sections %Q% (human60>0)) %Q% (strand == "+") %Q% (width > 99)))))
    ## next check end good MAPQ:
    if(length(mj))
        mj = as.character(seqnames(ends %Q% (seqnames %in% mj) %Q% (human60 > 0)))
    ## finally check seed qual:
    if(length(mj)){
        xs = x %Q% (seed) %Q% (rev(order(mapq))) %Q% (!duplicated(seqnames))
        cm = setNames(xs$mapq, xs$qname)
        if(any(!(mj %in% names(cm)))){
            cs = GRanges(paste0("chr", out.dt[qname %in% mj][!(qname %in% names(cm))]$peak), seqinfo=seqinfo(Hsapiens))
            cs$qname = out.dt[qname %in% mj][!(qname %in% names(cm))]$qname
            cs$seq = as.character(getSeq(Hsapiens, cs))
            ca = human[cs$seq] %Q% (rev(order(mapq))) %Q% (!duplicated(qname))
            cm = c(cm, setNames(as.integer(ca$mapq), cs[as.integer(ca$qname)]$qname))
        }
        mj = mj[cm[mj] == 60]
    }
    if(length(mj)){
        missedj = TRUE
        mate = reduce(parse.gr((x %Q% (seqnames %in% mj) %Q% (!seed) %Q% (mapq==60) %Q% (c_spec == "human"))$y) + 100)
        ## can't just be one -- in case of foldbacks break glass
        if(length(mate)>1){
            if(as.data.table(x %Q% (seqnames %in% mj) %Q% (!seed) %Q% (mapq==60) %Q% (c_spec == "human"))[width > 9][,.N,by=qname][,any(N>1)]){
                complex = TRUE
            } else if(!length((x %Q% (seqnames %in% mj) %Q% (seed) %Q% (mapq==60) %Q% (c_spec == "human")))){
                complex = TRUE
            } else{
                s = reduce(gr.flipstrand(parse.gr((x %Q% (seqnames %in% mj) %Q% (seed) %Q% (mapq==60) %Q% (c_spec == "human"))$y)))
                if(length(reduce(grbind(mate, s))) > 2){
                    complex = TRUE
                }
            }
        }
        inserts = gaps(sections %Q% (seqnames %in% mj) %Q% (human60>0 | seeds>0)) %Q% (strand == "+") %Q% (width > 19)
        if(length(inserts)){
            insertion = TRUE
            irep = x %&% inserts %Q% (c_spec!="human" | mapq < 60) %Q% (!seeds)
            if(any(irep$c_spec == "human" & irep$mapq < 60))
                irep[irep$c_spec == "human" & irep$mapq < 60]$c_spec = "low MAPQ"
            if(length(irep)){
                wides = as.data.table(gr.reduce(irep, by="query.id"))[, query.id[width > 19]]
                irep = irep %Q% (c_spec!="human") %Q% (c_type!="human" | hg19==score) %Q% (query.id %in% wides) %Q% (!duplicated(c_spec))
            }
        }
    }

    ## novel should not share section with seed!
    ## pulling out non-gaps of **human** alignments -- counting 10bp+ as a gap
    ## (ie not novel if this is true because the seed is bad)
    nh = setdiff(levels(seqnames(x)), unique(as.character(seqnames(gaps(sections %Q% (human>0)) %Q% (strand == "+") %Q% (width > 9)))))
    if(length(nh))
        nh = unique(as.character(seqnames(sections %Q% (seqnames %in% nh) %Q% (seeds == 0) %Q% (score > human60) %Q% (width > 9))))
    if(length(nh)){
        novel = TRUE
        nrep = c(nrep, x %&% (reduce(sections %Q% (seqnames %in% nh) %Q% (score > human60) %Q% (!seeds)) %Q% (width > 9)))
        if(length(nrep)){
            if(any(nrep$c_spec == "human" & nrep$mapq < 60))
                nrep[nrep$c_spec == "human" & nrep$mapq < 60]$c_spec = "low MAPQ"
            if(!"human60" %in% colnames(values(nrep))) nrep = nrep %$% sections
            if(insertion & any(!(nrep$human60))){
                nrep = nrep %Q% (human60)
            }
            if(length(nrep)){
                wides = as.data.table(gr.reduce(nrep, by="query.id"))[, query.id[width > 19]]
                nrep = nrep %Q% (c_spec!="human") %Q% (c_type!="human" | hg19==score) %Q% (query.id %in% wides) %Q% (!duplicated(c_spec))
            }
        }
    }

    ## no more that 19bp gap from *any* alignment
    nr = setdiff(levels(seqnames(x)), unique(as.character(seqnames(gaps(sections) %Q% (strand == "+") %Q% (width > 19)))))
    if(length(nr))
        nr = unique(as.character(seqnames(sections %Q% (seqnames %in% nr) %Q% (seeds == 0) %Q% (score > human60) %Q% (width > 19))))
    if(length(nr)){
        novel = TRUE
        nrep = c(nrep, x %Q% (seqnames %in% nr) %Q% (c_spec!="human" | mapq < 60) %Q% (!seeds))
        if(any(nrep$c_spec == "human"))
            nrep[nrep$c_spec == "human"]$c_spec = "low MAPQ"
        if(!"human60" %in% colnames(values(nrep))) nrep = nrep %$% sections
        if(insertion & any(!(nrep$human60))){
            nrep = nrep %Q% (human60)
        }
        if(length(nrep)){
            wides = as.data.table(gr.reduce(nrep, by="query.id"))[, query.id[width > 19]]
            nrep = nrep %Q% (c_spec!="human") %Q% (c_type!="human" | hg19==score) %Q% (query.id %in% wides) %Q% (!duplicated(c_spec))
        }
    }

    if(length(mj) == 0 & length(nh) == 0 & length(nr) == 0) {
        if(dt2gr(li) %N% uannot) reference = TRUE
        nrep = c(nrep, x %Q% (c_spec!="human" | mapq < 60) %Q% (!seeds))
        if(any(nrep$c_spec == "human" & nrep$mapq < 60))
            nrep[nrep$c_spec == "human" & nrep$mapq < 60]$c_spec = "low MAPQ"
        if(insertion & any(!(nrep$human60))){
            nrep = nrep %Q% (human60)
        }
        if(length(nrep)){
            wides = as.data.table(gr.reduce(nrep, by="query.id"))[, query.id[width > 19]]
            nrep = nrep %Q% (c_spec!="human") %Q% (c_type!="human" | hg19==score) %Q% (query.id %in% wides) %Q% (!duplicated(c_spec))
            novel = TRUE
        }
        g = gaps(sections) %Q% (strand == "+") %Q% (width == max(width))
        if(length(g)){
            if(width(g) / seqlengths(g)[as.character(seqnames(g))] > 0.75){
                g = g[width(g) / seqlengths(g)[as.character(seqnames(g))] > 0.75]
                mystery = TRUE
                una = x[1]
                una$c_spec = "unaligned sequence"
                nrep = c(nrep, una)
            } else{
                g = gaps(sections) %Q% (strand == "+") %Q% (width > 19)
                if(length(g)){
                    g = g[1]
                    una = x[1]
                    una$c_spec = "unaligned sequence"
                    nrep = c(nrep, una)
                    novel = TRUE
                }
            }
        }
    }

    if(missedj) refmap = novmap = TRUE
    mystery = mystery | !any(missedj, complex, reference, novel)
    return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap))
}

#' caller
#'
#' parses contig alignments to assign loose end to category and describe repeat types
#' @param li data.table loose end to evaluate
#' @param calns optional, data.table of contig alignments, default=NULL (will not parse alignments)
#' @param insert optional, integer pad representing insert length in bp to identify seed alignments, default=750
#' @param pad optional, window around loose end to allow contig seed windows, default=1e3
#' @param uannot optional, GRanges of unmappable annotations, default bin/101.unmappable.annotations.rds
#' @param ref_dir optional, path to directory of unzipped reference tarball, default assumes 'package/extdata/hg19_looseends'
#' @param ref_obj optional, list of BWA objects built from ref_dir fastas, names must match expected "human" "rep" "polyA" "microbe" (only "human" is used), default=NULL
#' @param return.contigs (logical)
#' @export
caller = function(li, calns = NULL, insert = 750, pad=NULL, uannot=NULL, ref_dir=system.file('extdata', 'hg19_looseends', package='loosends'), ref_obj=NULL, return.contigs = FALSE) {
    reference = FALSE
    complex = FALSE
    missedj = FALSE
    novel = FALSE
    mystery = FALSE
    insertion = FALSE
    refmap = FALSE
    novmap = FALSE
    rrep = nrep = irep = GRanges(c_spec=character())
    junction = ""
    if(is.null(pad)) pad = 1e3
    if(is.null(uannot)){
        uannot = readRDS(system.file('extdata', '101.unmappable.annotations.rds', package='loosends'))
    }
    if(!is.null(ref_obj)) { human = ref_obj$human
    } else{
        if(!file.exists(ref_dir)) stop("Provide correct ref_dir containing reference .fa files")
        if(!file.exists(paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37_decoy.fasta")
        human = BWA(fasta=paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))
    }
    
    if(is.null(calns) | !nrow(calns)) {
        nrep = uannot[1]
        nrep$c_spec = "mystery"
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            rrep = uannot %&% dt2gr(li)
            if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
        } else {
            refmap = TRUE
        }
        mystery = TRUE
        if (return.contigs) {
            res = list(call = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap),
                       irep = irep,
                       nrep = nrep,
                       contigs = GRanges())
            return(res)
        }
                       
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap))
    }
    if(inherits(calns$mapq, "character")) calns$mapq = as.integer(calns$mapq)
    calns[c_spec == "polyA", c_type := "rep"]
    calns[c_spec == "ribosomal", c_type := "ribosomal"]
    seed = GRanges(calns$peak)
    strand(seed) = ifelse(grepl("for", calns$track), "+", "-")
    seed$qname = calns$qname
    ## narrow down to 1kb window on correct strand
    ## this filters out some that were only caught on the reverse ....
    if("seed" %in% colnames(calns)) calns$seed = NULL
    calns = calns[!is.na(gr.match(seed, gr.flipstrand(dt2gr(li))+pad, ignore.strand=F))]
    if(!nrow(calns)){
        nrep = uannot[1]
        nrep$c_spec = "mystery"
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            rrep = uannot %&% dt2gr(li)
            if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
        } else {
            refmap = TRUE
        }
        mystery = TRUE
        if (return.contigs) {
            res = list(call = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap),
                       irep = irep,
                       nrep = nrep,
                       contigs = GRanges())
            return(res)
        }
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap))
    }
    
    ## filter out contigs with too few good input reads,
    ## unless there wouldn't be any left!
    ## keep going regardless!
    if(calns[map60.cov > 9, .N == 0]){
        reference = TRUE
        rrep = GRanges(li$seqnames, IRanges(li$start, li$end), c_spec="low MAPQ")
    }

    ## checking this isn't just the seed region fully mapping & not going anywhere
    ## also making sure it doesn't leave but come back
    seed = GRanges(calns$peak)
    strand(seed) = ifelse(grepl("for", calns$track), "+", "-")
    seed$qname = calns$qname
    calns$seed = !is.na(gr.match(dt2gr(calns), seed+insert, by="qname", ignore.strand=F))
    if(any(calns$seed) & all(calns[(seed), mapq == 60])) refmap = TRUE
    ch = cgChain(calns)
    x = ch@.galx
    values(x) = ch@values
    x$y = gr.string(dt2gr(calns)[floor(as.numeric(rownames(ch@values)))])
    ## ignoring indel cigar alignments -- w/r/t when width is used to prioritize repeats
    x$source.id = floor(as.numeric(rownames(ch@values)))
    ctl = si2gr(x)
    x = dt2gr(as.data.table(x)[, ":="(start = min(start), end = max(end)), by=source.id][!duplicated(source.id)][,!"source.id"])
    seqlengths(x) = seqlengths(ctl)[seqlevels(x)]
    
    ## now that we're in contig coordinates, parse for telomeres
    out.dt = calns[!duplicated(qname)]
    setkey(out.dt, qname)
    tgreps = find.telomeres(eighteenmer('g'), out.dt)
    tcreps = find.telomeres(eighteenmer('c'), out.dt)
    if(!is.null(tgreps)) x = merge.telomeres(tgreps, x, ch, out.dt, c_spec="G telomere")
    if(!is.null(tcreps)) x = merge.telomeres(tcreps, x, ch, out.dt, c_spec="C telomere")
    seqlengths(x) = seqlengths(ctl)[seqlevels(x)]
    ## viral calls overlapping telomere matches are ignored -- viral reference seqs include (TTAGGG)n
    ## actually viral calls overlapping anything human are ignored??
    if(any(x$c_type=="viral") & (!is.null(tgreps) | !is.null(tcreps))){
        virs = x %Q% (c_type == "viral")
        kvir = virs %O% (x %Q% (grepl("telomere", c_spec)))
        x = c(x %Q% (c_type != "viral"), virs[kvir < 0.9])
    }
    ## ribosomal calls overlapping polyA calls are ignored -- ribosomal ref includes (A)n
    if(any(x$c_spec=="ribosomal") & any((x %Q% (c_spec=="ribosomal")) %N% (x %Q% (c_spec=="polyA")))){
        ribs = x %Q% (c_spec == "ribosomal")
        krib = ribs %O% (x %Q% (c_spec == "polyA"))
        x = c(x %Q% (c_spec != "ribosomal"), ribs[krib < 0.9])
    }
    ## prioritze telomere motif matches over telomeric coordinate alignments
    if(any(x$c_spec=="Telomere / sub") & (!is.null(tgreps) | !is.null(tcreps))){
        tela = x %Q% (c_spec == "Telomere / sub")
        utel = is.na(gr.match(tela, x %Q% (grepl("tel", y)), by="seed"))
        x = c(x %Q% (c_spec != "Telomere / sub"), tela[utel])
    }

    
    ## make sure the contig goes and stays away:
    ## last 20 bases of the contig cannot align to seed window
    ends = gr.end(ctl, width=20)
    ends$seed = ends %N% (x %Q% (seed))
    ends$human60 = ends %N% (x %Q% (c_spec=="human") %Q% (mapq==60) %Q% (!(seed)))
    ## outp: qnames of contigs with a seed alignment in the last 20 bases
    outp = as.data.table(ends)[seed > 0, as.character(seqnames)]

    if(reference){
        rrep = copy(x[1])
        rrep$c_spec = "low MAPQ"
    }

    ## add to reference repeats: contig seed regions that also align elsewhere, at least 20bp
    ## (should be making sure the overlap specifically is at least 20bp)
    rrep = c(rrep, x %&% (x %Q% (seed)) %Q% (!(seed) | mapq<60 | c_type!="human") %Q% (width > 19) %Q% (seqnames %in% outp))
    if(any(rrep$c_spec == "human" & rrep$mapq < 60))
        rrep[rrep$c_spec == "human" & rrep$mapq < 60]$c_spec = "low MAPQ"
    rrep = rrep %Q% (c_spec!="human")
    if(length(rrep)) {
        reference = TRUE
    }

    ## filter out contigs with seed alignment in the last 20 bases
    x2 = x %Q% (!(seqnames %in% outp)) ## save x so you can check for reference repeats
    if(!length(x2)){
        nrep = uannot[1]
        nrep$c_spec = "mystery"
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            sup = uannot %&% dt2gr(li)
            if(is.na(sup$c_spec)) sup$c_spec = sup$repClass
            rrep = grbind(rrep, sup)
            rrep$c_spec = as.character(rrep$c_spec)
        } else refmap = refmap | all((x %Q% (seed))$mapq == 60)
        mystery = TRUE
        if (return.contigs) {
            res = list(call = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap),
                       irep = irep,
                       nrep = nrep,
                       contigs = x2)
            return(res)
        }
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap))
    } else x = x2
    
    ## trim the lowMQ bit to be only outside the repeat
##    if(any(x$c_type=="rep") & any((x %Q% (c_type=="rep")) %N% (x %Q% (c_spec=="human") %Q% (!seed & mapq < 60)))){
##        mq = x %Q% (!seed & c_spec == "human") %Q% (mapq < 60)
##        nmq = GenomicRanges::setdiff(mq, x %Q% (c_type=="rep"))
##        ov = gr.findoverlaps(nmq, mq, ignore.strand=F)
##        values(ov) = values(mq[ov$subject.id])
##        x = c(x %Q% (seed | c_spec != "human" | mapq==60), ov)
##    }

    ## remember you won't always have caught the junction -- might just be the other side
    sections = disjoin(x)
    sections$score = sections %N% x
    sections$seeds = sections %N% (x %Q% (seed))
    sections$seed60 = sections %N% (x %Q% (seed & mapq==60))
    sections$hg19 = sections %N% (x %Q% (c_type=="human"))
    sections$human = sections %N% (x %Q% (c_spec=="human"))
    sections$human60 = sections %N% (x %Q% (c_spec=="human") %Q% (mapq==60))

    x = x %**% sections
    if(any(x$c_type=="viral")){
        virs = x %Q% (c_type == "viral")
        vwid = as.data.table(virs)[, sum(width), by=.(query.id, seed, qname)]
        nh = (x %Q% (c_type != "human" | mapq < 60)) %Q% (c_type != "viral")
        nwid = as.data.table(nh)[, sum(width), by=.(query.id, seed,qname)][, max(V1), by=.(seed, qname)]
        nwid[, sq := paste(seed, qname)]
        setkey(nwid, sq)
        vwid[, nwid := nwid[.(paste(vwid$seed, vwid$qname)), V1]]
        vwid[is.na(nwid), nwid := 0]
        virs = virs %Q% (query.id %in% vwid[V1 > nwid][, query.id])
        x = c(x %Q% (c_type != "viral"), virs)

        sections = disjoin(x)
        sections$score = sections %N% x
        sections$seeds = sections %N% (x %Q% (seed))
        sections$seed60 = sections %N% (x %Q% (seed & mapq==60))
        sections$hg19 = sections %N% (x %Q% (c_type=="human"))
        sections$human = sections %N% (x %Q% (c_spec=="human"))
        sections$human60 = sections %N% (x %Q% (c_spec=="human") %Q% (mapq==60))
    }
    
    x = as.data.table(x)
    x[, clength := seqlengths(ctl)[qname]]

    if(x[, any(seed)] & x[(seed), all(mapq<60)]){
        reference = TRUE
        rrep = rrep %Q% (!(seqnames %in% outp))
        if(length(setdiff(colnames(values(dt2gr(x))), colnames(values(rrep)))))
            rrep = rrep %$% dt2gr(x)[, setdiff(colnames(values(dt2gr(x))), colnames(values(rrep)))]
        rrep = grbind(rrep, dt2gr(copy(x[(seed)][!duplicated(qname)])[, c_spec := "low MAPQ"], seqlengths=seqlengths(ctl)))
        rrep$c_spec = as.character(rrep$c_spec)
    }
    x = dt2gr(x)
    seqlengths(x) = seqlengths(ctl)[seqlevels(x)]
    sections = sections %$% (x[, 'clength']) %$% (x %Q% (c_spec != "human"))[, 'c_spec']
    sections = as.data.table(sections)
    
    sections[, seed.end := ifelse(any(seeds>0), max(end[seeds>0], na.rm=T), as.integer(0)), by=seqnames]
    sections[, seed.side := end <= seed.end]
    sections[, mate.side := start > seed.end]

    sections[, seed.map := ifelse(seed.end<19, "mystery", ifelse(any(seed60 == seeds & seed60 > 0 & seed.side & width>19), "mappable", "unmappable")), by=seqnames]
    sections[is.na(seed.map), seed.map := "mystery"]
    sections[, mate.map := ifelse(seed.end+19 >= clength, "mystery", ifelse(sum(width[human60 > 0 & mate.side]) >= clength-seed.end-19, "mappable", "unmappable")), by=seqnames]
    sections[is.na(mate.map), mate.map := "mystery"]
    if(length(ends %Q% (seed == 0) %Q% (human60 > 0))){
        me = as.character(seqnames(ends %Q% (seed == 0) %Q% (human60 > 0)))
        sections[seqnames %in% me, mate.map := ifelse(sum(width[human60 > 0 & mate.side]) >= clength-seed.end-99, "mappable", "unmappable"), by=seqnames]
    }

    sm = sections[, any(seed.map!="mystery" & mate.map!="mystery"), by=seqnames][(V1), as.character(seqnames)]
    if(length(sm)){
        sections = sections[seqnames %in% sm]
        x = x %Q% (seqnames %in% sm)
        seqlevels(x) = sm
    }

    cs = sections[, sum(width[seed.side & !seeds]) > 99, by=seqnames][(V1), as.character(seqnames)]
    if(length(cs)){
        if(all(sections[, seqnames %in% cs])){
            mystery = TRUE
            nrep = dt2gr(sections[(mate.side)][1][, c_spec := "mystery"])
            novmap = FALSE
            if (return.contigs) {
                res = list(call = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap),
                           irep = irep,
                           nrep = nrep,
                           contigs = x)
                return(res)
            }
            return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap, junction=junction))
        }
        sections = sections[!(seqnames %in% cs)]
        x = x %Q% (!(seqnames %in% cs))
        seqlevels(x) = setdiff(seqlevels(x), cs)
    }
    
    if(sections[,all(seed.map=="mystery")]){
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            sup = (uannot %&% dt2gr(li) %Q% (width == max(width)))[1]
            rrep = x %Q% (!duplicated(qname))
            rrep$c_spec = sup$c_spec
            rrep$c_type = "rep"
            if(is.na(sup$c_spec)) rrep$c_spec = sup$repClass
            sections[, seed.map := "unmappable"]
        } else{
            cs = out.dt[qname %in% as.character(sections$seqnames), GRanges(paste0("chr", peak), seqinfo=seqinfo(Hsapiens), qname=qname)]
            cs$seq = as.character(getSeq(Hsapiens, cs))
            ca = human[cs$seq] %Q% (rev(order(mapq))) %Q% (!duplicated(qname))
            cm = setNames(as.integer(ca$mapq), cs[as.integer(ca$qname)]$qname)
            sections[, seed.map := ifelse(cm[as.character(seqnames)]==60, "mappable", "unmappable")]
            if(any(cm < 60)){
                sup = out.dt[qname %in% names(cm)[cm<60]]
                sup$c_type = "human"
                sup$c_spec = "low MAPQ"
                rrep = c(rrep, dt2gr(sup, seqlengths = seqlengths(x)))
            }
        }
    }    

    sections = dt2gr(as.data.table(sections))
    seqlengths(sections) = seqlengths(ctl)[seqlevels(sections)]
    

    ## can only catch a reference repeat this way if the seed is part of the contig
    rr = unique(as.character(seqnames(x %Q% (seeds>0) %Q% (score > 1 | human>human60) %Q% (width > 19))))
    if(length(rr)){
        reference = TRUE
        rrep = rrep %Q% (seqnames %in% as.character(seqnames(x)))
        if(length(setdiff(colnames(values(x)), colnames(values(rrep)))))
            rrep = rrep %$% x[, setdiff(colnames(values(x)), colnames(values(rrep)))]
        rrep = c(rrep, x %Q% (seqnames %in% rr) %Q% (seeds>0) %Q% (!(seed) | seeds > 1 | human>human60) %Q% (width > 19))
        rrep = rrep %Q% (!(c_spec=="human" & mapq==60))
        if(any(rrep$c_spec == "human" & rrep$mapq < 60))
            rrep[rrep$c_spec == "human" & rrep$mapq < 60]$c_spec = "low MAPQ"
    }    

    ## make sure if there are multiple they go to the same place..
    ## seed has good MAPQ -- if not included in any contig, align using getSeq
    ## last 20 bases overlap good MAPQ
    ## no MAPQ60 gaps 100bp+
    ## first check gaps:
    gapcs = gaps(sections) %Q% (strand=="+") %Q% (width > 19)
    sections$gap = as.character(seqnames(sections)) %in% as.character(seqnames(gapcs))
    if(any(sections$gap & sections$mate.map=="mappable")) sections[sections$gap & sections$mate.map=="mappable"]$mate.map = "unmappable"
    
    mj = setdiff(levels(seqnames(x)), unique(as.character(seqnames(gaps(sections %Q% (human60>0)) %Q% (strand == "+") %Q% (width > 99)))))
    mj = as.data.table(sections)[seqnames %in% mj, any(seed.map=="mappable" & mate.map=="mappable"), by=seqnames][(V1), as.character(seqnames)]
    ## next check end good MAPQ:
    if(length(mj))
        mj = as.character(seqnames(ends %Q% (seqnames %in% mj) %Q% (human60 > 0)))
    ## finally check seed qual:
    if(length(mj)){
        xs = x %Q% (seed) %Q% (rev(order(mapq))) %Q% (!duplicated(seqnames))
        cm = setNames(xs$mapq, xs$qname)
        if(length(cm) & all(cm < 60)) mj = mj[0]
        if(any(!(mj %in% names(cm)))){
            cs = GRanges(paste0("chr", out.dt[qname %in% mj][!(qname %in% names(cm))]$peak), seqinfo=seqinfo(Hsapiens))
            cs$qname = out.dt[qname %in% mj][!(qname %in% names(cm))]$qname
            cs$seq = as.character(getSeq(Hsapiens, cs))
            ca = human[cs$seq] %Q% (rev(order(mapq))) %Q% (!duplicated(qname))
            cm = c(cm, setNames(as.integer(ca$mapq), cs[as.integer(ca$qname)]$qname))
        }
        mj = mj[cm[mj] == 60]
    }
    if(length(mj)){
        missedj = TRUE
        refmap = TRUE
        novmap = TRUE

        rrep = rrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% mj)
        reference = length(rrep) > 0
        nrep = nrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% mj)
        novel = length(nrep) > 0
        sections = sections %Q% (seqnames %in% mj)
        seqlevels(sections) = mj
        x = x %Q% (seqnames %in% mj)
        seqlevels(x) = mj

        mate = reduce(parse.gr((x %Q% (seqnames %in% mj) %Q% (!seed) %Q% (mapq==60) %Q% (c_spec == "human"))$y) + 100)
        junction = paste(c(gr.string(dt2gr(li)), gr.string(gr.start(mate, ignore.strand=F))), collapse = " | ")
        ## can't just be one -- in case of foldbacks break glass
        if(length(mate)>1){
            if(as.data.table(x %Q% (seqnames %in% mj) %Q% (!seed) %Q% (mapq==60) %Q% (c_spec == "human"))[width > 9][,.N,by=qname][,any(N>1)]){
                complex = TRUE
            } else if(!length((x %Q% (seqnames %in% mj) %Q% (seed) %Q% (mapq==60) %Q% (c_spec == "human")))){
                complex = TRUE
            } else{
                s = reduce(gr.flipstrand(parse.gr((x %Q% (seqnames %in% mj) %Q% (seed) %Q% (mapq==60) %Q% (c_spec == "human"))$y)))
                if(length(reduce(grbind(mate, s))) > 2){
                    complex = TRUE
                }
            }
        }
        inserts = gaps(sections %Q% (seqnames %in% mj) %Q% (human60>0 | seeds>0)) %Q% (strand == "+") %Q% (width > 19)
        if(length(inserts)){
            insertion = TRUE
            irep = x %&% inserts %Q% (c_spec!="human" | mapq < 60) %Q% (!seeds)
            if(any(irep$c_spec == "human" & irep$mapq < 60))
                irep[irep$c_spec == "human" & irep$mapq < 60]$c_spec = "low MAPQ"
##            if(length(irep)){
##                wides = as.data.table(gr.reduce(irep, by="query.id"))[, query.id[width > 19]]
##                irep = irep %Q% (c_spec!="human") %Q% (c_type!="human" | hg19==score) %Q% (query.id %in% wides) %Q% (!duplicated(c_spec))
##            }
        }
    }

    rrep = c(rrep, x %Q% (seeds>0) %Q% (c_spec!="human" | mapq < 60))
    nrep = c(nrep, x %Q% (seeds==0) %Q% (c_spec!="human" | mapq < 60))
    if(any(nrep$c_spec == "human" & nrep$mapq < 60))
        nrep[nrep$c_spec == "human" & nrep$mapq < 60]$c_spec = "low MAPQ"
    if(length(nrep)){
        if(any(is.na(nrep$query.id))){
            max.qi = max(nrep$query.id, na.rm=T)
            if(!is.finite(max.qi)) max.qi = 0
            nrep[is.na(nrep$query.id)]$query.id = 1:sum(is.na(nrep$query.id)) + max.qi
        }
        si = seqinfo(nrep)
        nrep = gr.reduce(nrep, by="query.id") %Q% (width > 19)
        seqinfo(nrep) = si[seqlevels(nrep)]
        if(length(nrep)){
            novel = TRUE
        }
    }

    if(any(rrep$c_spec == "human" & rrep$mapq < 60))
        rrep[rrep$c_spec == "human" & rrep$mapq < 60]$c_spec = "low MAPQ"
    if(length(rrep)){
        if(any(is.na(rrep$query.id))){
            max.qi = max(rrep$query.id, na.rm=T)
            if(!is.finite(max.qi)) max.qi = 0
            rrep[is.na(rrep$query.id)]$query.id = 1:sum(is.na(rrep$query.id)) + max.qi
        }
        si = seqinfo(rrep)
        rrep = gr.reduce(rrep, by="query.id") %Q% (width > 19)
        seqinfo(rrep) = si[seqlevels(rrep)]
        if(length(rrep)){
            reference = TRUE
        }
    }

    nr = as.character( seqnames(sections %Q% (seed.map == "mappable"))) %>% unique
    if(length(nr)){
        refmap = TRUE
        rrep = rrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% nr)
        reference = length(rrep) > 0
        nrep = nrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% nr)
        novel = length(nrep) > 0

        sections = sections %Q% (seqnames %in% nr)
        seqlevels(sections) = nr
        x = x %Q% (seqnames %in% nr)
        seqlevels(x) = nr
    }

    nh = as.character(seqnames(sections %Q% (mate.map == "mappable"))) %>% unique
    if(length(nh)){
        novmap = TRUE
        rel = sections %Q% (seqnames %in% nh)
        if(all(rel$seed.map == "mystery")){
            if(dt2gr(li) %N% uannot){
                reference = TRUE
                rrep = uannot %&% dt2gr(li)
                if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
            } else {
                refmap = TRUE
            }
        } else{
            refmap = any(rel$seed.map=="mappable")
        }
        nrep = nrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% nh)
        novel = length(nrep) > 0
        if(any(nrep$c_spec == "human" & nrep$mapq < 60))
            nrep[nrep$c_spec == "human" & nrep$mapq < 60]$c_spec = "low MAPQ"
        if(!"human60" %in% colnames(values(nrep))) nrep = nrep %$% sections        
        if(insertion & any(!(nrep$human60))){
            nrep = nrep %Q% (human60)
        }

        sections = sections %Q% (seqnames %in% nh)
        seqlevels(sections) = nh
        x = x %Q% (seqnames %in% nh)
        seqlevels(x) = nh
    }

    g = gaps(sections) %Q% (strand == "+") %Q% (width > 19)
    if(length(g)){
        g = g[1]
        una = x[1]
        una$c_spec = "unaligned sequence"
        nrep = c(nrep, una)
        novel = TRUE
    }

    if(missedj) refmap = novmap = TRUE
    mystery = mystery | !any(missedj, complex, reference, novel)

    ## browser()

    call.dt = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap, junction=junction)

    if (return.contigs) {
        res = list(contigs = x, call = call.dt, irep = irep, nrep = nrep)
        return(res)
    }
    return(call.dt)
    ##return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap, junction=junction))
}


#' gr.sum.strand
#'
#' strand-specific implementing gr.sum
gr.sum.strand = function(gr){
    pos = gr.sum(gr %Q% (strand == "+"))
    strand(pos) = "+"
    neg = gr.sum(gr %Q% (strand == "-"))
    strand(neg) = "-"
    return(c(pos, neg))
}

#' read.based
#'
#' classifying loose end based on discordant read alignments by constructing pseudo-contigs based on consensus alignment patterns
#' @param li data.table loose end to classify
#' @param ri data.table read alignments
#' @param pad optional, padding around loose end li to search for discordant reads, default=1e3
#' @param ref_dir optional, path to directory of unzipped reference tarball, default assumes 'package/extdata/hg19_looseends'
#' @param ref_obj optional, list of BWA objects built from ref_dir fastas, names must match expected "human" "rep" "polyA" "microbe" (only "human" is used), default=NULL
#' @param return.full also return reads
#' @param 
read.based = function(li, ri, pad=NULL, ref_dir=system.file('extdata', 'hg19_looseends', package='loosends'), ref_obj=NULL){
    if(is.null(pad)) pad = 1e3
    t = ifelse(li$strand == "+", "sample.rev", "sample.for")
    ri = ri[rev(order(qwidth))][!duplicated(paste(R1, qname))]
    pp = (gr.tile(gr.flipstrand(dt2gr(li))+pad, 200) %Q% (width == 200))[,c()]

    out.dt = rbindlist(lapply(1:length(pp), function(i){
        win = pp[i]
        ri[, seed := dt2gr(ri) %N% win > 0 & track==t]        
        rtmp = ri[qname %in% qname[seed]]
        rtmp[(seed)][is.dup(qname), seed := R1 == R1[1], by=qname] ## in case of a foldback where both reads are "seed", choose one as the anchor
        rtmp[is.dup(paste(qname, R1)), strand := ifelse(seed | !any(seed), strand, ifelse(strand=="-", "+", "-"))] ## in case of split alignment OF SEED READ, flip the strand of non-seed split (keep split of mate)
        rtmp = rtmp[mapq == 60]
        if(nrow(rtmp)==0 | rtmp[(seed), .N==0] | rtmp[!(seed), .N==0]) return(NULL)
        ctig = reduce(gr.sum.strand(dt2gr(rtmp[(seed)])) %Q% (score > 0))
        ctig = ctig[ctig %NN% dt2gr(rtmp[(seed)]) > 9]
        if(length(ctig)==0) return(NULL)
        mate = reduce(gr.sum.strand(gr.flipstrand(dt2gr(rtmp[!(seed)]))) %Q% (score > 0))
        mate = mate[mate %NN% gr.flipstrand(dt2gr(rtmp[!(seed)])) > 9]
        seed = mate %NN% ctig
        ctig = as.data.table(reduce(c(gr.fix(ctig, mate), gr.fix(mate[seed > 0], ctig))))
        ctig[, read.cov := dt2gr(ctig) %NN% dt2gr(rtmp[(seed)])] ##rtmp[, sum(seed)]]
        ctig = rbind(ctig, as.data.table(mate[seed==0])[, read.cov := rtmp[, sum(!seed)]], use.names=T, fill=T)
        ctig[, map60.cov := read.cov]
        ctig[, marker := c(0, cumsum(width))[1:.N]]
        ctig[, me := sum(width) - cumsum(width)]
        ctig[, cigar := paste0(ifelse(marker > 0, paste0(marker, "S"), ""), width, "M", ifelse(me > 0, paste0(me, "S"), ""))]
        ctig[strand=="-", cigar := paste0(ifelse(me > 0, paste0(me, "S"), ""), width, "M", ifelse(marker > 0, paste0(marker, "S"), ""))]
        ctig[, peak := gr.string(gr.stripstrand(win[,c()]))]
        ctig[, qname := i]
        return(ctig[, !c("marker", "me")])
    }), use.names=TRUE, fill=TRUE)

    if(nrow(out.dt)){
        out.dt[, track := t]
        out.dt[, leix := li$leix]
        out.dt[, qname := paste(li$sample, leix, track, qname, sep="_")]
        out.dt[, c_type := "human"]
        out.dt[, c_spec := "human"]
        out.dt$mapq = 60
        out.dt$seq = ""
    } else{
        out.dt$qname = character()
        out.dt$track = character()
        out.dt$leix = character()
        out.dt$c_type = character()
        out.dt$c_spec = character()
        out.dt$mapq = integer()
        out.dt$seq = character()
    }

    dt = caller(li, out.dt[is.dup(qname)], ref_dir=ref_dir, ref_obj=ref_obj)
    return(dt[, c("complex", "missedj", "junction")])
}

#' transform
#'
#' this function allows conversion between data.table and gr.sum
#' (used for convenience with a data.table by= argument)
#' @param seq character seqnames
#' @param s integer start
#' @param e integer end
transform = function(seq, s, e){
    gr = GRanges(seq, IRanges(s, e))
    dt = gr2dt(gr.sum(gr))
    dt$seqnames = as.character(dt$seqnames)
    return(dt)
}

#' ggraph.loose.ends
#'
#' this function processes a single gGraph to identify
#' high quality loose ends and categorize them
#' @param gg a gGraph built from a JaBbA model (ie with loose ends)
#' @param cov.rds character path to binned coverage file (should contain GRanges)
#' @param tbam character path to tumor BAM file
#' @param nbam optional, character path to normal BAM file, default=NULL
#' @param id optional, character sample id, default=NULL
#' @param outdir optional, character path to output directory, default=NULL (will not write output files if outdir=NULL)
#' @param purity optional, fractional purity of sample, default assumes 1
#' @param ploidy optional, ploidy of sample, default infers from gg
#' @param field optional, name of informative column in cov.rds, default="ratio"
#' @param PTHRESH optional, threshold for GLM p-value for calling true positive loose ends, default=3.4e-7 provides consanguinity with large dataset bonferroni correction
#' @param verbose optional, will print status messages, default=FALSE
#' @param mc.cores optional, parallel cores for building contigs per loose end, default=1
#' @param ref_dir optional, path to directory of unzipped reference tarball, default assumes 'package/extdata/hg19_looseends'
#' @param overwrite, optional, logical indicating whether to generate new output files if corresponding files have already been written to outdir, default=FALSE (will load existing files)
#' @export
ggraph.loose.ends = function(gg, cov.rds, tbam, nbam=NULL, id=NULL, outdir=NULL, purity=NULL, ploidy=NULL, field="ratio", PTHRESH=3.4e-7, verbose=F, mc.cores=1, ref_dir=system.file('extdata', 'hg19_looseends', package='loosends'), ref_obj=NULL, overwrite=FALSE){
    if(!is.null(outdir)) if(!file.exists(outdir)) {
                             tryCatch(readLines(pipe(paste("mkdir", outdir))), error = function(e) stop(paste("Provided output directory", outdir, "does not exist and cannot be made")))
                         }
    if(is.null(ref_obj)){
        if(!file.exists(ref_dir)) stop("Provide correct ref_dir containing reference .fa files")
        if(!file.exists(paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37_decoy.fasta")
        if(!file.exists(paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))) stop("ref_dir must contain msilab_combined_TraFicv8-3_satellites.fa")
        if(!file.exists(paste(ref_dir, "PolyA.fa", sep="/"))) stop("ref_dir must contain PolyA.fa")
        if(!file.exists(paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37.withviral.fasta")
    }

    le.all = filter.graph(gg, cov.rds, purity=purity, ploidy=ploidy, field=field, PTHRESH=PTHRESH, verbose=verbose)
    process.loose.ends(le.all[(true.pos)], tbam, nbam=nbam, id=id, outdir=outdir, mc.cores=mc.cores, ref_dir=ref_dir, ref_obj=ref_obj, verbose=verbose, overwrite=overwrite)
}

#' process.loose.ends
#'
#' this function categorizes an input set of loose ends
#' runs process.single.end for each input loose end
#' @param le data.table or GRanges of loose ends to analyze
#' @param tbam character path to tumor BAM file
#' @param nbam optional, character path to normal BAM file, default=NULL
#' @param id optional, character sample id, default=NULL
#' @param outdir optional, character path to output directory, default=NULL (will not write output files if outdir=NULL)
#' @param mc.cores optional, parallel cores for building contigs per loose end, default=1
#' @param ref_dir optional, path to directory of unzipped reference tarball, default assumes 'package/extdata/hg19_looseends'
#' @param verbose optional, default=FALSE
#' @param overwrite, optional, logical indicating whether to generate new output files if corresponding files have already been written to outdir, default=FALSE (will load existing files)
#' @export
process.loose.ends = function(le, tbam, nbam=NULL, id=NULL, outdir=NULL, mc.cores=1, ref_dir=system.file('extdata', 'hg19_looseends', package='loosends'), ref_obj=NULL, verbose=FALSE, overwrite=FALSE){
    le = as.data.table(copy(le))
    if(is.null(ref_obj)){
        if(!file.exists(ref_dir)) stop("Provide correct ref_dir containing reference .fa files")
        if(!file.exists(paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37_decoy.fasta")
        if(!file.exists(paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))) stop("ref_dir must contain msilab_combined_TraFicv8-3_satellites.fa")
        if(!file.exists(paste(ref_dir, "PolyA.fa", sep="/"))) stop("ref_dir must contain PolyA.fa")
        if(!file.exists(paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37.withviral.fasta")
        if(!(length(tbam) == length(nbam) & length(nbam) == length(id))) stop("tbam, nbam, and id must all be length=1 or length=length(le)")
        if(length(tbam) > 1 & length(tbam) != nrow(le)) stop("if tbam, nbam, and id are length>1, must be length(le)")
        human = BWA(fasta=paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))
        rep = BWA(fasta=paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))
        polyA = BWA(fasta=paste(ref_dir, "PolyA.fa", sep="/"))
        microbe = BWA(fasta=paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"), keep_sec_with_frac_of_primary_score=0.2)
        ref_obj = list(human=human, rep=rep, polyA=polyA, microbe=microbe)
    }
    if(length(tbam)==1){
        tbam = rep(tbam, nrow(le))
        nbam = rep(nbam, nrow(le))
        id = rep(id, nrow(le))
    }
    rbindlist(mclapply(1:nrow(le), function(i) process.single.end(le[i], tbam[i], nbam=nbam[i], id=id[i], outdir=outdir, ref_dir=ref_dir, ref_obj=ref_obj, verbose=verbose, overwrite=overwrite)[, i := i], mc.cores=mc.cores), fill=T, use.names=T)
}

#' process.single.end
#'
#' this function loads reads surrounding a loose end and their mates,
#' assembles reads into contigs, aligns contigs to references,
#' and parses alignments to categorize the loose end
#' @param li GRanges or data.table equivalent representing loose end to evaluate
#' @param tbam character string representing path to tumor BAM file (must have corresponding BAM index file)
#' @param nbam optional, character string representing path to normal BAM file (if used, must have corresponding BAM index file)
#' @param id optional, character string representing sample identifier
#' @param outdir optional, character string representing path to output directory
#' if provided, four output files are created
#'         leix.reads.rds -- contains data.table of realigned read pairs originating within 5 kbp of loose end
#'         leix.contigs.rds -- contains data.table of assembled contigs, one row per contig
#'         leix.aligned.contigs.rds -- contains data.table of contig alignments
#'         leix.call.rds -- contains one row data.table describing categorization of provided loose end li (same as returned value)
#' @param ref_dir optional, path to directory of unzipped reference tarball, default assumes 'package/extdata/hg19_looseends'
#' @param ref_obj optional, list of BWA objects built from ref_dir fastas, names must match expected "human" "rep" "polyA" "microbe", default=NULL
#' @param verbose optional, default=FALSE
#' @param overwrite, optional, logical indicating whether to generate new output files if corresponding files have already been written to outdir, default=FALSE (will load existing files)
#' @return one row data table describing categorization of provided loose end li
#' @export
process.single.end = function(li, tbam, nbam=NULL, id=NULL, outdir=NULL, ref_dir=system.file('extdata', 'hg19_looseends', package='loosends'), ref_obj=NULL, verbose=FALSE, overwrite=FALSE){
    li = as.data.table(copy(li))
    message(gr.string(dt2gr(li)))
    if(nrow(li) > 1) stop("More than 1 loose end provided; did you mean 'process.loose.ends'?")
    if(nrow(li)==0) stop("0 loose ends provided")
    if(!is.null(outdir)) if(!file.exists(outdir)) {
                             tryCatch(readLines(pipe(paste("mkdir", outdir))), error = function(e) stop(paste("Provided output directory", outdir, "does not exist and cannot be made")))
                         }
    if(is.null(ref_obj)){
        if(!file.exists(ref_dir)) stop("Provide correct ref_dir containing reference .fa files")
        if(!file.exists(paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37_decoy.fasta")
        if(!file.exists(paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))) stop("ref_dir must contain msilab_combined_TraFicv8-3_satellites.fa")
        if(!file.exists(paste(ref_dir, "PolyA.fa", sep="/"))) stop("ref_dir must contain PolyA.fa")
        if(!file.exists(paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37.withviral.fasta")

        human = BWA(fasta=paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))
        rep = BWA(fasta=paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))
        polyA = BWA(fasta=paste(ref_dir, "PolyA.fa", sep="/"))
        microbe = BWA(fasta=paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"), keep_sec_with_frac_of_primary_score=0.2)
        ref_obj = list(human=human, rep=rep, polyA=polyA, microbe=microbe)
    } else{
        human = ref_obj$human
        rep = ref_obj$rep
        polyA = ref_obj$polyA
        microbe = ref_obj$microbe
    }
    uannot = readRDS(system.file('extdata', '101.unmappable.annotations.rds', package='loosends'))
    
    if(is.null(id)) id = "SAMPLE"
    if(!"sample" %in% colnames(li)) li$sample = id
    if(!"leix" %in% colnames(li)) li$leix = paste(id, strsplit(gr.string(dt2gr(li)), "-")[[1]][1], sep=":")
    ro = !is.null(outdir)

    if(ro & file.exists(paste(outdir, paste0(li$leix, ".call.rds"), sep="/"))){
        call = readRDS(paste(outdir, paste0(li$leix, ".call.rds"), sep="/"))
    } else{

        leix = li$leix
        if(ro & file.exists(paste(outdir, paste0(li$leix, ".reads.rds"), sep="/"))) {
            ri = readRDS(paste(outdir, paste0(li$leix, ".reads.rds"), sep="/"))
        } else{
            
            ri = loose.reads(li, tbam=tbam, nbam=nbam, filter=FALSE, pad=5e3, ref=human, verbose=verbose)
            ri$sample = as.character(ri$sample)
            ri$leix = li$leix
            ri[, track := paste(ifelse(sample == li$sample, "sample", "control"), ifelse(strand == "+", "for", "rev"), sep=".")]
            ri[, concord := !(loose.pair) & .N == 2 & length(unique(seqnames)) == 1 & strand[R1] != strand[R2] & strand[start == min(start)]=="+" & min(start) + 3e3 > max(start), by=qname]
            ri[, anchor := (loose.pair & high.mate) | ( !(loose.pair) & mapq > 50 & !(concord))]
            ri$seq = as.character(ri$seq)
            ri$start = as.integer(ri$start)
            ri$end = as.integer(ri$end)
            ri$flag = as.integer(ri$flag)
            if(!("reading.frame" %in% colnames(ri))){
                ri$reading.frame = ri$seq
                ri[strand == "-", reading.frame := as.character(reverseComplement(DNAStringSet(reading.frame)))]
            }

            if(ro) saveRDS(ri, paste(outdir, paste0(li$leix, ".reads.rds"), sep="/"))
            gc()
        }

        somatic = as.logical(nrow(ri[grepl("control", track)]))
        wholseed = dt2gr(li)[,c()]+1e3

        pp = gr.stripstrand(gr.tile(wholseed, 200) %Q% (width == 200))

        .build.tigs = function(ri, pp, ro, outdir){
            if(verbose) message("assembling contigs")
            out.dt = rbindlist(lapply(1:length(pp), function(i){
                win = pp[i]
                build.from.win(win, ri)
            }), use.names=TRUE, fill=TRUE)
            if("seed" %in% colnames(ri)) ri$seed = NULL
            out.dt$sample = rep(id, nrow(out.dt))
            out.dt$leix = rep(leix, nrow(out.dt))
            out.dt$somatic = rep(somatic, nrow(out.dt))
            out.dt = out.dt[!is.na(seq)]

            if(verbose) message(paste("assembled", nrow(out.dt), "contigs across all tracks"))
            if(nrow(out.dt)){
                out.dt[, qname := 1:.N, by=.(track, sample, leix)]
                out.dt[, qname := paste(sample, leix, track, qname, sep="_")]
            } else{
                out.dt$qname = character()
            }
            if(verbose) message("aligning contigs to 4 references")
            ureps = rep[out.dt$seq]
            preps = polyA[out.dt$seq]
            mreps = microbe[out.dt$seq]
            hreps = human[out.dt$seq]
            virals = seqlevels(microbe)[85:length(seqlevels(microbe))]
            
            if(!(length(ureps))){
                ureps$cigar = character()
                ureps$mapq = character()
                ureps$AS = integer()
                ureps$flag = character()
            }
            if(!(length(preps))){
                preps$cigar = character()
                preps$mapq = character()
                preps$AS = integer()
                preps$flag = character()
            }
            if(!(length(hreps))){
                hreps$cigar = character()
                hreps$mapq = character()
                hreps$AS = integer()
                hreps$flag = character()
            }
            if(!(length(mreps))){
                mreps$cigar = character()
                mreps$mapq = character()
                mreps$AS = integer()
                mreps$flag = character()
            }
            keep.cols = c("cigar", "flag", "mapq", "AS")
            values = rbind(out.dt[as.integer(ureps$qname)], out.dt[as.integer(preps$qname)], out.dt[as.integer(hreps$qname)], fill=T, use.names=T)
            ralns = rbind(as.data.table(ureps[, keep.cols]), as.data.table(preps[, keep.cols]), as.data.table(hreps[, keep.cols]))
            ralns = cbind(ralns, values)
            if(nrow(out.dt)){
                rch = cgChain(ralns)
                good.ids = c(as.character(seqnames(gaps(rch$x) %Q% (strand=="+"))), setdiff(out.dt$qname, ralns$qname))
                mreps$query.id  = as.integer(mreps$qname)
                mreps$qname = out.dt[mreps$query.id, qname]
                vreps = mreps %Q% (seqnames %in% virals) %Q% (qname %in% good.ids)
                if(verbose & length(vreps)) message("adding viral alignments")
                valns = cbind(as.data.table(vreps[, keep.cols]), out.dt[vreps$query.id])
                calns = rbind(ralns, valns)
                calns[, c_type := c(rep("rep", length(ureps)), rep("polyA", length(preps)), rep("human", length(hreps)), rep("viral", length(vreps)))]
                calns[, c_spec := c_type]
                calns[c_type == "rep", c_spec := dunlist(strsplit(as.character(seqnames), "#", fixed=T))[rev(!duplicated(rev(listid)))][order(as.integer(listid)), V1]]
            } else{
                calns = ralns
            }
            calns$somatic = rep(somatic, nrow(calns))
            calns$mapq = as.integer(calns$mapq)
            if(nrow(out.dt)){
                if(verbose) message("quantifying contig read support")
                if(nrow(valns)){
                    ch = cgChain(calns)
                } else ch = rch
                rs = getSeq(Hsapiens, trim(gr.fix(gr.chr(wholseed + 15000), Hsapiens)))
                refseq = as.character(rs)
                supp.dat = rbindlist(lapply(1:nrow(out.dt), function(i){
                    op = FALSE
                    if(grepl("control", out.dt[i]$track)) return(list(support=0, cov=0, mapq60=0, used.cs=op))
                    ctigs = calns[qname == out.dt[i]$qname]
                    if(nrow(ctigs)==0) return(list(support=0, cov=0, mapq60=0, used.cs=op))
                    win = out.dt[i, GRanges(peak)]
                    seed = ri[dt2gr(ri) %N% win > 0]
                    strand(win) = ifelse(grepl("for", out.dt[i]$track), "+", "-")
                    rtmp = ri[qname %in% seed$qname]
                    x = dt2gr(as.data.table(ch$x)[seqnames == out.dt[i]$qname])
                    if(any(x == reduce(x)) | nrow(ctigs) == 1){ 
                        rc = BWA(seq=c(ctigs[1]$seq, refseq))[rtmp$reading.frame] %Q% (seqnames == 1) %Q% (mapq == 60)
                        rc$qname = rtmp[as.integer(rc$qname), qname]
                    } else{
                        ## change this to legacy contig support for now
                        rc = contig.support.legacy(dt2gr(rtmp), ctigs, refseq)
                        op=TRUE
                    }
                    if(!length(rc)) return(list(support=0, cov=0, mapq60=0, used.cs=op))
                    t = rtmp[qname %in% rc$qname, table(factor(grepl("sample", track), c(TRUE, FALSE)))]
                    supp = ri[qname %in% rc$qname][grepl("sample", track)]
                    return(data.table(support=t["TRUE"] / sum(t), cov=seed[grepl("sample", track) & qname %in% rc$qname, .N], mapq60=seed[grepl("sample", track) & qname %in% rc$qname & mapq==60, .N], used.cs=op))
                }), use.names=T)
                out.dt$support = supp.dat$support
            } else {
                out.dt$support = numeric()
                out.dt$cov = numeric()
                out.dt$mapq60 = numeric()
            }
            setkey(out.dt, qname)
            calns[, tumor.spec := qname %in% out.dt[support==1]$qname]
            calns$support = out.dt[.(calns$qname), support]
            calns$read.cov = out.dt[.(calns$qname), cov]
            calns$map60.cov = out.dt[.(calns$qname), mapq60]

            if(ro) saveRDS(out.dt, paste(outdir, paste0(li$leix, ".contigs.rds"), sep="/"))
            if(verbose) message("parsing contigs for telomeric matches")
            all.contigs = copy(calns)
            all.contigs[c_type == "human", c_spec := ifelse(seqnames %in% c(1:22, "X", "Y"), "human", "unassembled")]
            all.contigs$cmer = munch(all.contigs, eighteenmer('c'))
            all.contigs$gmer = munch(all.contigs, eighteenmer('g'))

            if(verbose) message("parsing alignments for unmappable repeat overlaps")
            coord.calls = copy(all.contigs[c_type=="human"])[, !c("c_type", "c_spec")]
            ov = gr.findoverlaps(dt2gr(coord.calls), uannot)
            if(length(ov)){
                subj = ov$subject.id
                strand(ov) = coord.calls[ov$query.id, strand]
                values(ov) = values(dt2gr(coord.calls)[ov$query.id])
                ov$c_type = "rep"
                ov$c_spec = uannot[subj]$c_spec
                coord.calls = as.data.table(ov)
                if(ro) saveRDS(rbind(all.contigs, coord.calls, fill=T, use.names=T), paste(outdir, paste0(li$leix, ".aligned.contigs.rds"), sep="/"))
                return(rbind(all.contigs, coord.calls, fill=T, use.names=T))
            }
            if(ro) saveRDS(all.contigs, paste(outdir, paste0(li$leix, ".aligned.contigs.rds"), sep="/"))
            return(all.contigs)
        }

        if(ro & file.exists(paste(outdir, paste0(li$leix, ".aligned.contigs.rds"), sep="/"))){
            all.contigs = readRDS(paste(outdir, paste0(li$leix, ".aligned.contigs.rds"), sep="/"))
        } else{
            all.contigs = .build.tigs(ri, pp, ro, outdir)
        }

        if(verbose) message("generating call...")
        call = caller(li, all.contigs, ref_obj=ref_obj)
        disc = read.based(li, ri, ref_obj=ref_obj)
        recall = (!call$missedj & disc$missedj) | (!call$complex & disc$complex)
        call[, missedj := missedj | disc$missedj]
        call[, complex := complex | disc$complex]
        call[junction == "" & disc$junction != "", junction := disc$junction]
        call[junction != "" & disc$junction != "", junction := {
            cgr = parse.gr(junction)
            dgr = parse.gr(disc$junction)
            s = gr.reduce(cgr[1], dgr[1], ignore.strand=F)
            m = gr.reduce(cgr[-1], dgr[-1], ignore.strand=F)
            m = gr.start(reduce(m + 200) - 200, ignore.strand=F)
            paste(c(gr.string(s), gr.string(m)), collapse=" | ")
        }]
        call[, mystery := !missedj & !complex & mate.mappable & seed.mappable & !insertion]
        if(any(recall)){
            call[recall]$call = update.call(call[recall])
            call[(missedj), seed.mappable := TRUE]
            call[(missedj), mate.mappable := TRUE]
        }
        if(call$mystery){
            if(verbose) message("mystery: repeating assembly at larger intervals")
            pp = gr.stripstrand((gr.tile(wholseed-250, 500)+250) %Q% (width == 1e3))
            wide.contigs = .build.tigs(ri, pp, ro, outdir)
            call = caller(li, wide.contigs, ref_obj=ref_obj)
        }
        if(ro) saveRDS(call, paste(outdir, paste0(li$leix, ".call.rds"), sep="/"))
    }
    
    call[, category := ifelse(complex, "complex rearrangement", ifelse(missedj, "missed junction", ifelse(mystery | grepl("mystery", mate.repeats), "mystery", paste("type", as.integer(!seed.mappable) + as.integer(!mate.mappable), "loose end"))))]
    call[category=="mystery", mystery := TRUE]
    call[, ":="(
        leix = li[, leix],
        loose.end = li[, paste0(seqnames, ":", start, strand)],
        sample = id)]
    gc()
    return(call)
}

#' .sample.spec
#'
#' loads reads and mates for a single sample (tumor or normal)
#' @param le GRanges or data.table of loose ends
#' @param bam path to BAM file
#' @param pad integer width of padding to add around loose ends
#' @param verbose optional, default=FALSE
.sample.spec = function(le, bam, pad, verbose=FALSE){
    if(verbose) message(paste("loading reads from", bam))
    if(!inherits(le, "GRanges")) le = dt2gr(le)
    sl = c(seqlengths(BamFile(bam)), setNames(1, "*"))
    has.chr = any(grepl("chr", seqnames(seqinfo(BamFile(bam)))))
    if(has.chr) { w = gr.chr(le) + pad
    } else w = le + pad
    reads = as.data.table(unlist(bamUtils::read.bam(bam, gUtils::gr.reduce(gr.stripstrand(w)), pairs.grl=T, isDuplicate=NA, isPaired=TRUE, tag="SA")))
    splits = reads[!is.na(SA)]
    if(nrow(splits) > 0){
        splits$SA = as.character(splits$SA)
        splwin = dunlist(strsplit(splits$SA, ";"))
        spl = unlist(lapply(strsplit(splwin$V1, ","), function(w) paste(w[1], w[2], sep=":")))
        spl = GRanges(spl)
        spl$qname = splits[as.integer(splwin$listid)]$qname
        splitsides = as.data.table(unlist(read.bam(bam, gUtils::gr.reduce(spl+150)-150, pairs.grl=T, isDuplicate=NA, tag="SA")) %Q% (qname %in% spl$qname))[order(mrnm, mpos)][!duplicated(paste(seqnames, start, qname, seq))]
        reads = rbind(reads, splitsides, fill=T, use.names=TRUE)
    }
    reads[, unpmate := bamflag(flag)[, "hasUnmappedMate"]==1]
    reads[, isunp := start == 1 & is.na(seq)]
    reads[, unp := any(unpmate) & any(isunp), by=qname]
    reads[(unp), ":="(start = ifelse(isunp, start[unpmate], start), end = ifelse(isunp, end[unpmate], end)), by=qname]
    reads[, missing := any(is.na(seq)), by=qname]
    mw = reads[!is.na(mrnm) & !is.na(seq)]
    stopifnot(!is.na(mw$mrnm))
    mw[, ":="(seqnames = mrnm, start = ifelse(is.na(mpos), start, mpos), end = ifelse(is.na(mpos), end, mpos))]
    reads = reads[!is.na(seq)]
    ##    fixmr = reads[qname %in% mw$qname, setNames(mrnm, qname)]
    ##    mw[, seqnames := fixmr[qname]]
    mw = dt2gr(mw[,!"strand"], seqlengths=sl)
    mw[width(mw) == 0] = mw[width(mw) == 0] + 1
    mate.wins = gUtils::gr.reduce(mw+150)-150
    reads = reads[, !c("unpmate", "isunp", "unp", "SA")]
    if(length(mate.wins) > 0){
        if(verbose) message(paste("loading mates from", length(mate.wins), "windows"))
        m.mate.wins = mate.wins        
        mates = rbindlist(lapply(seq(1, length(m.mate.wins), 100), function(i){
            mi = read.bam(bam, gUtils::gr.reduce(gr.stripstrand(m.mate.wins[i:min(length(m.mate.wins), i+99)])), pairs.grl=F, isDuplicate=NA)
            mi = as.data.table(mi[mi$qname %in% reads$qname])
            gc()
            return(mi)
        }), fill=TRUE, use.names=TRUE)
        gc()
        if(verbose) message("mates loaded")
        rpair = rbind(reads, mates, fill=TRUE, use.names=TRUE)[!duplicated(paste(qname, flag)),]; rpair$MQ = NULL
    } else {
        rpair = reads[!duplicated(paste(qname, flag)),]; rpair$MQ = NULL
    }
    rpair[, R1 := bamflag(flag)[, "isFirstMateRead"]==1]
    rpair[, R2 := bamflag(flag)[, "isSecondMateRead"]==1]
    rpair[, paired := any(R1) & any(R2), by=qname]
    if(verbose) message(ifelse(all(rpair$paired), "Found All Mates!!", "Some mates still missing - perhaps BAM was deduplicated"))
    rpair[, MQ := rev(mapq), by=qname]
    rpair[, count := .N, by = qname]
    rpair[count == 0, MQ := 0]
    ## rpair[, MQ := MQ * as.integer(.N>1), by=qname]
    ##    flip = rpair$strand == "-"
    flip = bamflag(rpair$flag)[, "isMinusStrand"] == 1 | rpair$strand == "-"
    ##flip = bamflag(rpair$flag)[, "isMinusStrand"] == 1
    rpair[flip, seq := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq)))]
    ##    rpair = rpair[!duplicated(paste(qname, R1, seq))]
    rpair = rpair[rev(order(nchar(seq)))][!duplicated(paste(qname, R1))]
    reads = dt2gr(rpair)
    gc()
    return(reads)
}

#' .realign
#'
#' realigns reads and their mates single-end to attach individual MAPQs
#' @param reads GRanges or data.table of reads from BAM file with $seq in original reading frame
#' @param ref reference genome for realignment
#' @param gg optional, gGraph of sample used to identify sequences fitted in graph, default=NULL
#' @param filter optional, filter=TRUE will return loose read pairs only, filter=FALSE returns all, default=TRUE
#' @param verbose optional, default=FALSE
.realign = function(reads, ref, gg=NULL, filter=TRUE, verbose=F){
    if(!is.null(gg)){
        seqs = unique(seqnames(gg$nodes[!is.na(cn)]$gr))
    } else {
        seqs = c(1:22, "X", "Y")
        if(any(grepl("chr", seqlevels(ref)))) seqs = gr.chr(seqs)
    }
    if(filter){
        qni = as.data.table(reads)[is.na(mapq) | mapq<50 | is.na(MQ), (unique(qname))]
    } else{
        qni = unique(reads$qname)
    }
    redo = reads$qname %in% qni
    if(verbose) message(paste("realigning", sum(redo), "reads"))
    realn = ref[setNames(reads$seq, 1:length(reads))[redo]]
    realn$mapq = as.integer(realn$mapq)
    realn$query.id = as.integer(realn$qname); realn$qname = reads[realn$query.id]$qname
    values(realn) = cbind(values(realn), values(reads[realn$query.id, !(colnames(values(reads)) %in% colnames(values(realn)))]))
    uix = ! ( (1:length(reads))[redo] %in% realn$query.id )
    uix = (1:length(reads))[redo][uix]
    unaln = as.data.table(reads[uix])
    unaln[, ":="(
        seqnames = "*",
        start = 1,
        end = 0,
        flag = ifelse(R1, 69, 113),
        mapq = 0,
        query.id = uix)]
    realn = rbind(as.data.table(realn), unaln, fill=T, use.names=TRUE)
    realn$reading.frame = reads[realn$query.id]$seq
    realn$mapq = as.integer(realn$mapq)
    realn[, mapq := ifelse(!(seqnames %in% seqs), as.integer(0), mapq), by=query.id]
    realn = realn[rev(order(mapq))][!duplicated(query.id), ]
    gc()
    realn[, MQ := ifelse(rep(.N, .N)==1, as.integer(NA), c(mapq[-1], mapq[1])), by=qname]
    cols = c(colnames(realn)[colnames(realn) %in% colnames(as.data.table(reads))], "reading.frame", "AS")
    realn = realn[, cols, with=F]
    lqn = realn[mapq>50 & (is.na(MQ) | MQ < 1), qname]
    if(filter){
        realn = realn[qname %in% lqn,]
        realn[, loose.pair := TRUE]
    } else realn[, loose.pair := qname %in% lqn]
    realn[, high.mate := mapq>50 & (is.na(MQ) | MQ < 1)]
    gc()
    return(realn)
}

#' loose.reads
#'
#' loads and reads and their mates from given windows and realigns single end for read-specific MAPQ
#' @param le GRanges or data.table windows around which to sesarch for reads
#' @param tbam character path to tumor BAM file
#' @param pad optional, integer padding around le windows, default=25e3
#' @param nbam optional, character path to normal BAM file, default=NULL
#' @param ref optional, BWA object of reference genome for realignment or path to fasta, default=package/extdata/hg19_looseends/human_g1k_v37_decoy.fasta
#' @param filter optional, logical filter=TRUE returns loose read pairs only, filter=FALSE returns all read pairs annotated with logical $loose column, default=T
#' @param gg optional, gGraph corresponding to sample, used to identify sequences fitted in graph, default=NULL
#' @param verbose optional, default=FALSE
#' @export
loose.reads = function(le, tbam, pad=25e3, nbam=NA, ref=system.file('extdata', 'hg19_looseends', 'human_g1k_v37_decoy.fasta', package='loosends'), filter=TRUE, gg=NULL, verbose=FALSE){
    le = copy(le)
    if(inherits(ref, "character")){
        if(!file.exists(ref)) stop("Provide reference BWA object or path to reference fasta for loose.reads")
        ref = RSeqLib::BWA(fasta=ref)
    }
    if(is.na(nbam)) nbam = NULL
    id = le[1]$sample
    treads = .sample.spec(copy(le), tbam, pad, verbose=verbose)
    realn = .realign(treads, ref, filter=filter, gg=gg, verbose=verbose)
    realn$sample = id

    if(!is.null(nbam)){
        nreads = .sample.spec(copy(le), nbam, pad, verbose=verbose)
        nrealn = .realign(nreads, ref, filter=filter, gg=gg, verbose=verbose)
        nrealn$sample = paste0(id, "N")
        gc()
        return(rbind(realn, nrealn, fill=TRUE, use.names=TRUE))
    }
    gc()
    return(realn)
}

#' @name prep_loose_reads
#' @title prep_loose_reads
#'
#' Get loose reads data table ready for assembly
#'
#' @param li
#' @param loose.reads.dt
#'
#' @return data.table
#' @export
prep_loose_reads = function(li, loose.reads.dt) {
    if (nrow(li) > 1) {
        stop("More than 1 loose end provided; did you mean 'process.loose.ends'?")
    }
    if (nrow(li)==0) {
        stop("0 loose ends provided")
    }

    if (is.null(li$sample)) {
        stop("li must contain column $sample")
    }

    if (is.null(li$leix)) {
        stop("li must contain column $leix")
    }

    if (is.null(loose.reads.dt$sample)) {
        stop("loose.reads.dt must contain column $sample")
    }

    ## get loose ends objects ready
    ri = copy(loose.reads.dt)
    ri$sample = as.character(ri$sample)

    ## add leix to both ri and li because this is needed for saving samples
    ri$leix = li$leix

    ## denote sample vs. control
    ri[, track := paste(ifelse(sample == li$sample, "sample", "control"), ifelse(strand == "+", "for", "rev"), sep=".")]
    ri[, concord := !(loose.pair) & .N == 2 & length(unique(seqnames)) == 1 & strand[R1] != strand[R2] & strand[start == min(start)]=="+" & min(start) + 3e3 > max(start), by=qname]
    ri[, anchor := (loose.pair & high.mate) | ( !(loose.pair) & mapq > 50 & !(concord))]
    ri$seq = as.character(ri$seq)
    ri$start = as.integer(ri$start)
    ri$end = as.integer(ri$end)
    ri$flag = as.integer(ri$flag)
    if(!("reading.frame" %in% colnames(ri))){
        ri$reading.frame = ri$seq
        ri[strand == "-", reading.frame := as.character(reverseComplement(DNAStringSet(reading.frame)))]
    }

    return(ri)
}

#' @name assemble_loose_reads
#' @title assemble_loose_reads
#'
#' this function loads reads surrounding a loose end and their mates,
#' assembles reads into contigs, aligns contigs to references,
#' and parses alignments to categorize the loose end
#' @param loose.end.gr GRanges
#' @param loose.reads.dt data.table output from loose.reads
#' @param verbose optional, default=FALSE
#' @return data.table with assembled contigs
#' @export
assemble_loose_reads = function(li,
                                loose.reads.dt,
                                ref_dir=system.file('extdata', 'hg19_looseends', package='loosends'),
                                ref_obj=NULL,
                                verbose=FALSE,
                                overwrite=FALSE) {
    li = as.data.table(copy(li))
    if (verbose) {
        message("Loose end coordinates: ", gr.string(dt2gr(li)))
    }

    if(is.null(ref_obj)){
        if(!file.exists(ref_dir)) stop("Provide correct ref_dir containing reference .fa files")
        if(!file.exists(paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37_decoy.fasta")
        if(!file.exists(paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))) stop("ref_dir must contain msilab_combined_TraFicv8-3_satellites.fa")
        if(!file.exists(paste(ref_dir, "PolyA.fa", sep="/"))) stop("ref_dir must contain PolyA.fa")
        if(!file.exists(paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37.withviral.fasta")

        human = BWA(fasta=paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))
        rep = BWA(fasta=paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))
        polyA = BWA(fasta=paste(ref_dir, "PolyA.fa", sep="/"))
        microbe = BWA(fasta=paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"), keep_sec_with_frac_of_primary_score=0.2)
        ref_obj = list(human=human, rep=rep, polyA=polyA, microbe=microbe)
    } else{
        human = ref_obj$human
        rep = ref_obj$rep
        polyA = ref_obj$polyA
        microbe = ref_obj$microbe
    }
    uannot = readRDS(system.file('extdata', '101.unmappable.annotations.rds', package='loosends'))


    ri = copy(loose.reads.dt)

    somatic = as.logical(nrow(ri[grepl("control", track)]))
    wholseed = dt2gr(li)[,c()]+1e3

    pp = gr.stripstrand(gr.tile(wholseed, 200) %Q% (width == 200))


    message(li$leix)
    message(li$id)
    all.contigs = .build.tigs(ri, pp, li$sample, li$leix, FALSE, NULL)
    return(all.contigs)
}

#' @name call_loose_end
#' @title call_loose_end
#'
#' @description
#'
#' This function calls a loose end given reads from near to the loose end
#'
#' @param li (GRanges) loose end ranges, must have metatdata $leix and $sample
#' @param ri (data.table) loose reads table prepped (e.g. from prep_loose_reads)
#' @param ref_dir (character) path to directory with fasta/bwa indices
#' @param ref_obj (character) or a list of BWA indices, if ref_dir not provided
#' @param pad (numeric) window size for local assembly, default 1 kbp
#' @param mix.tn (logical) mix tumor/normal reads before assembly
#' @param verbose
#'
#' @return list with entries
#' - all.contigs (data.table)
#' - wide.contigs (data.table or NULL)
#' - call (data.table)
call_loose_end = function(li, ri,
                          ref_dir = NULL,
                          ref_obj = NULL,
                          pad = 1e3,
                          mix.tn = FALSE,
                          verbose = FALSE) {
    
    if(is.null(ref_obj)){
        if(!file.exists(ref_dir)) stop("Provide correct ref_dir containing reference .fa files")
        if(!file.exists(paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37_decoy.fasta")
        if(!file.exists(paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))) stop("ref_dir must contain msilab_combined_TraFicv8-3_satellites.fa")
        if(!file.exists(paste(ref_dir, "PolyA.fa", sep="/"))) stop("ref_dir must contain PolyA.fa")
        if(!file.exists(paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37.withviral.fasta")

        human = BWA(fasta=paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))
        rep = BWA(fasta=paste(ref_dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))
        polyA = BWA(fasta=paste(ref_dir, "PolyA.fa", sep="/"))
        microbe = BWA(fasta=paste(ref_dir, "human_g1k_v37.withviral.fasta", sep="/"), keep_sec_with_frac_of_primary_score=0.2)
        ref_obj = list(human=human, rep=rep, polyA=polyA, microbe=microbe)
    } else{
        human = ref_obj$human
        rep = ref_obj$rep
        polyA = ref_obj$polyA
        microbe = ref_obj$microbe
    }
    
    uannot = readRDS(system.file('extdata', '101.unmappable.annotations.rds', package='loosends'))

    id = li$sample


    if (verbose) {
        message("Generating call")
    }

    .build.tigs = function(ri, pp, id, leix, verbose = TRUE) {
        if(verbose) message("assembling contigs")
        out.dt = rbindlist(lapply(1:length(pp), function(i){
            win = pp[i]
            build.from.win(win, ri, verbose = verbose)
        }), use.names=TRUE, fill=TRUE)
        if("seed" %in% colnames(ri)) ri$seed = NULL
        out.dt$sample = rep(id, nrow(out.dt))
        out.dt$leix = rep(leix, nrow(out.dt))
        out.dt$somatic = rep(somatic, nrow(out.dt))
        out.dt = out.dt[!is.na(seq)]

        if(verbose) message(paste("assembled", nrow(out.dt), "contigs across all tracks"))
        if(nrow(out.dt)){
            out.dt[, qname := 1:.N, by=.(track, sample, leix)]
            out.dt[, qname := paste(sample, leix, track, qname, sep="_")]
        } else{
            out.dt$qname = character()
        }
        if(verbose) message("aligning contigs to 4 references")
        ureps = rep[out.dt$seq]
        preps = polyA[out.dt$seq]
        mreps = microbe[out.dt$seq]
        hreps = human[out.dt$seq]
        virals = seqlevels(microbe)[85:length(seqlevels(microbe))]
        
        if(!(length(ureps))){
            ureps$cigar = character()
            ureps$mapq = character()
            ureps$AS = integer()
            ureps$flag = character()
        }
        if(!(length(preps))){
            preps$cigar = character()
            preps$mapq = character()
            preps$AS = integer()
            preps$flag = character()
        }
        if(!(length(hreps))){
            hreps$cigar = character()
            hreps$mapq = character()
            hreps$AS = integer()
            hreps$flag = character()
        }
        if(!(length(mreps))){
            mreps$cigar = character()
            mreps$mapq = character()
            mreps$AS = integer()
            mreps$flag = character()
        }
        keep.cols = c("cigar", "flag", "mapq", "AS")
        values = rbind(out.dt[as.integer(ureps$qname)], out.dt[as.integer(preps$qname)], out.dt[as.integer(hreps$qname)], fill=T, use.names=T)
        ralns = rbind(as.data.table(ureps[, keep.cols]), as.data.table(preps[, keep.cols]), as.data.table(hreps[, keep.cols]))
        ralns = cbind(ralns, values)
        if(nrow(out.dt)){
            rch = cgChain(ralns)
            good.ids = c(as.character(seqnames(gaps(rch$x) %Q% (strand=="+"))), setdiff(out.dt$qname, ralns$qname))
            mreps$query.id  = as.integer(mreps$qname)
            mreps$qname = out.dt[mreps$query.id, qname]
            vreps = mreps %Q% (seqnames %in% virals) %Q% (qname %in% good.ids)
            if(verbose & length(vreps)) message("adding viral alignments")
            valns = cbind(as.data.table(vreps[, keep.cols]), out.dt[vreps$query.id])
            calns = rbind(ralns, valns)
            calns[, c_type := c(rep("rep", length(ureps)), rep("polyA", length(preps)), rep("human", length(hreps)), rep("viral", length(vreps)))]
            calns[, c_spec := c_type]
            calns[c_type == "rep", c_spec := dunlist(strsplit(as.character(seqnames), "#", fixed=T))[rev(!duplicated(rev(listid)))][order(as.integer(listid)), V1]]
        } else{
            calns = ralns
        }
        calns$somatic = rep(somatic, nrow(calns))
        calns$mapq = as.integer(calns$mapq)
        if(nrow(out.dt)){
            if(verbose) message("quantifying contig read support")
            if(nrow(valns)){
                ch = cgChain(calns)
            } else ch = rch
            rs = getSeq(Hsapiens, trim(gr.fix(gr.chr(wholseed + 15000), Hsapiens)))
            refseq = as.character(rs)
            supp.dat = rbindlist(lapply(1:nrow(out.dt), function(i){
                op = FALSE
                if(grepl("control", out.dt[i]$track)) return(list(support=0, cov=0, mapq60=0, used.cs=op))
                ctigs = calns[qname == out.dt[i]$qname]
                if(nrow(ctigs)==0) return(list(support=0, cov=0, mapq60=0, used.cs=op))
                win = out.dt[i, GRanges(peak)] ## should this window be padded?
                seed = ri[dt2gr(ri) %N% win > 0] ## which qnames overlap the peak window?
                strand(win) = ifelse(grepl("for", out.dt[i]$track), "+", "-")
                rtmp = ri[qname %in% seed$qname] ## which reads correspond with that qname?
                x = dt2gr(as.data.table(ch$x)[seqnames == out.dt[i]$qname])
                if(any(x == reduce(x)) | nrow(ctigs) == 1){ 
                    rc = BWA(seq=c(ctigs[1]$seq, refseq))[rtmp$reading.frame] %Q% (seqnames == 1) %Q% (mapq == 60)
                    rc$qname = rtmp[as.integer(rc$qname), qname]
                } else{
                    ## change to legacy contig support
                    rc = contig.support.legacy(dt2gr(rtmp), ctigs, refseq)
                    op=TRUE
                }
                if(!length(rc)) {
                    return(list(support=0, cov=0, mapq60=0, used.cs=op))
                }
                ## reannotate reads with track information
                ## t = rtmp[qname %in% rc$qname, table(factor(grepl("sample", track), c(TRUE, FALSE)))]
                ## supp = ri[qname %in% rc$qname][grepl("sample", track)]
                ## return(data.table(support=t["TRUE"] / sum(t),
                ##                   cov=seed[grepl("sample", track) & qname %in% rc$qname, .N],
                ##                   mapq60=seed[grepl("sample", track) & qname %in% rc$qname & mapq==60, .N],
                ##                   used.cs=op))
                ## how many reads are from the tumor?
                n.sample.reads = rtmp[sample == id & qname %in% rc$qname, .N]
                n.supp.reads = rtmp[qname %in% rc$qname, .N]
                return(data.table(support= n.sample.reads / n.supp.reads,
                                  cov=seed[grepl("sample", track) & qname %in% rc$qname, .N],
                                  mapq60=seed[grepl("sample", track) & qname %in% rc$qname & mapq==60, .N],
                                  used.cs=op))                                                    
            }), use.names=T)
            out.dt$support = supp.dat$support
        } else {
            out.dt$support = numeric()
            out.dt$cov = numeric()
            out.dt$mapq60 = numeric()
        }
        setkey(out.dt, qname)
        calns[, tumor.spec := qname %in% out.dt[support==1]$qname]

        ## also keep the proportion of tumor-specific supporting reads
        calns$support = out.dt[.(calns$qname), support]

        ## get seed coverage and mapq60 from this
        calns$read.cov = out.dt[.(calns$qname), cov]
        calns$map60.cov = out.dt[.(calns$qname), mapq60]

        if(verbose) message("parsing contigs for telomeric matches")
        all.contigs = copy(calns)
        if (nrow(all.contigs)) {
            all.contigs[c_type == "human", c_spec := ifelse(seqnames %in% c(1:22, "X", "Y"), "human", "unassembled")]
            all.contigs$cmer = munch(all.contigs, eighteenmer('c'))
            all.contigs$gmer = munch(all.contigs, eighteenmer('g'))

            if(verbose) message("parsing alignments for unmappable repeat overlaps")
            coord.calls = copy(all.contigs[c_type=="human"])[, !c("c_type", "c_spec")]
            ov = gr.findoverlaps(dt2gr(coord.calls), uannot)
            if(length(ov)){
                subj = ov$subject.id
                strand(ov) = coord.calls[ov$query.id, strand]
                values(ov) = values(dt2gr(coord.calls)[ov$query.id])
                ov$c_type = "rep"
                ov$c_spec = uannot[subj]$c_spec
                coord.calls = as.data.table(ov)
                return(rbind(all.contigs, coord.calls, fill=T, use.names=T))
            }
        }
        return(all.contigs)
    }

    if (verbose) {
        message("Building contigs")
    }
    somatic = as.logical(nrow(ri[grepl("control", track)]))
    wholseed = dt2gr(li)[,c()] + pad ##dt2gr(li)[,c()]+1e3
    if (verbose) {
        message("Using window size: ", pad, " bp")
    }
    pp = gr.stripstrand(gr.tile(wholseed, 200) %Q% (width == 200))

    ## check if contigs should be assembled mixing tumor and normal reads
    ri.input = copy(ri)
    if (mix.tn) {
        if (verbose) {
            message("Mixing tumor/normal reads before assembly")
        }
        ri.input$track = NULL ## set track to null so that assembly is sample-agnostic
    }
    all.contigs = .build.tigs(ri.input, pp, li$sample, li$leix, verbose = TRUE)

    ## also save the filtered contigs
    if (verbose) {
        message("Generating call")
    }

    if (mix.tn) {
        res = caller(li, all.contigs, ref_obj = ref_obj, return.contigs = TRUE)
    } else {
        ## if not mixed, call should only be generated from tumor contigs
        res = caller(li, all.contigs[grepl('sample', track),], ref_obj = ref_obj, return.contigs = TRUE)
    }
    
    call = res$call
    filtered.contigs = res$contigs
    disc = read.based(li, ri, ref_obj=ref_obj)
    recall = (!call$missedj & disc$missedj) | (!call$complex & disc$complex)
    call[, missedj := missedj | disc$missedj]
    call[, complex := complex | disc$complex]
    call[junction == "" & disc$junction != "", junction := disc$junction]
    call[junction != "" & disc$junction != "", junction := {
        cgr = parse.gr(junction)
        dgr = parse.gr(disc$junction)
        s = gr.reduce(cgr[1], dgr[1], ignore.strand=F)
        m = gr.reduce(cgr[-1], dgr[-1], ignore.strand=F)
        m = gr.start(reduce(m + 200) - 200, ignore.strand=F)
        paste(c(gr.string(s), gr.string(m)), collapse=" | ")
    }]
    call[, mystery := !missedj & !complex & mate.mappable & seed.mappable & !insertion]
    if(any(recall)){
        call[recall]$call = update.call(call[recall])
        call[(missedj), seed.mappable := TRUE]
        call[(missedj), mate.mappable := TRUE]
    }
    if(call$mystery){
        if(verbose) message("mystery: repeating assembly at larger intervals")
        pp = gr.stripstrand((gr.tile(wholseed-250, 500)+250) %Q% (width == 1e3))
        ## changed from ri to ri.input for consistency
        wide.contigs = .build.tigs(ri.input, pp, li$sample, li$leix, verbose = TRUE)
        if (mix.tn) {
            call = caller(li, wide.contigs, ref_obj=ref_obj)
        } else {
            call = caller(li, wide.contigs[grepl('sample', track),], ref_obj = ref_obj)
        }
    }

    call[, category := ifelse(complex, "complex rearrangement", ifelse(missedj, "missed junction", ifelse(mystery | grepl("mystery", mate.repeats), "mystery", paste("type", as.integer(!seed.mappable) + as.integer(!mate.mappable), "loose end"))))]
    call[category=="mystery", mystery := TRUE]
    call[,  ":="(
        leix = li[, leix],
        loose.end = li[, paste0(seqnames, ":", start, strand)],
        sample = id)]

    if (verbose) {
        message("category: ", call$category)
    }

    res = list(all.contigs = all.contigs, call = call, wide.contigs = data.table(), filtered.contigs = as.data.table(filtered.contigs))
    if (exists('wide.contigs')) {
        if (inherits(wide.contigs, "data.table")) {
            res$wide.contigs = wide.contigs
        }
    }
    return(res)
}

#' @name read_support
#' @title read_support
#'
#' @description
#'
#' wrapper around contig.support for easy(-ish) debugging
#'
#' like contig.support, the input contigs must be nonempty and contain a single qname
#'
#' @param loose.end (GRanges)
#' @param loose.end.str (character)
#' @param reads.dt (data.table) data.table of reads coercible to GRanges
#' @param contigs.dt (data.table) data.table with colnames track, peak, qname. all qnames should be identical.
#' @param pad (numeric) pad loose ends for getting the refseq?
#' @param seed.pad (numeric) pad around
#' @param all.reads (logical) return all input reads or just the chunks supporing the contig? (default FALSE)
#' @param ref.coordinates (logical) want reads in ref (vs contig) coordinates? default TRUE
#' @param verbose (logical)
#' @param ... additional parameters to contig.support
#' 
#' @return data.table of candidate reads
#' if all.reads is FALSE and there are no supporting reads an empty table is returned
#' if all.reads is TRUE the table will have column $supporting indicating whether read suupports given contig
read_support = function(loose.end = NA, loose.end.str = NA_character_,
                        reads.dt = data.table(),
                        contigs.dt = data.table(),
                        ref.bwa = NULL,
                        fasta = NA_character_,
                        pad = 1e3,
                        refseq.pad = 1e4,
                        seed.pad = 0,
                        chimeric = TRUE,
                        strict = TRUE,
                        min.bases = 20,
                        min.aligned.frac = 0.95,
                        isize.diff = 1e3,
                        all.reads = FALSE,
                        ref.coordinates = TRUE,
                        verbose = FALSE,
                        ...) {

    if (!is.na(loose.end) & inherits(loose.end, 'GRanges')) {
        li = copy(loose.end)
    } else if (!is.na(loose.end.str)) {
        li = parse.gr(loose.end.str)
    } else {
        stop("valid loose end not provided")
    }

    if (!length(li)) {
        stop("valid loose end not provided")
    }

    if (!nrow(reads.dt)) {
        stop("Empty reads.dt")
    }
    
    if (!nrow(contigs.dt)) {
        stop("Either contigs.dt or junction must be supplied")
    }

    if (length(unique(contigs.dt$qname)) > 1) {
        stop("some contigs have different qnames")
    }

    if (is.null(ref.bwa)) {

        if (file.exists(fasta)) {
            if (verbose) {
                message("Building BWA object from supplied FASTA: ", fasta)
            }
            ref.bwa = BWA(fasta = fasta)
        } else {
            if (verbose) {
                message("Grabbing reference sequence around loose end")
            }
            wholseed = li[,c()] + pad
            rs = getSeq(Hsapiens, trim(gr.fix(gr.chr(wholseed + 15000), Hsapiens)))
            refseq = as.character(rs)
            ref.bwa = BWA(refseq)
        }
        
    } else {
        if (inherits(ref.bwa, 'BWA')) {
            if (verbose) {
                message("Using supplied ref.bwa")
            }
        } else {
            stop("Invalid format supplied for BWA")
        }
    }

    if (verbose) {
        message("Grabbing reads overlapping contig peak")
    }
    win = contigs.dt[1, GRanges(peak)]
    seed = reads.dt[dt2gr(reads.dt) %N% (win + seed.pad) > 0] ## which qnames overlap the peak window?
    strand(win) = ifelse(grepl("for", contigs.dt[1]$track), "+", "-")
    window.reads.dt = reads.dt[qname %in% seed$qname] ## which reads correspond with that qname?

    if (verbose) {
        message("Starting contig support")
    }
    rc = contig.support(dt2gr(window.reads.dt),
                        contigs.dt,
                        ref = ref.bwa,## refseq,
                        chimeric = chimeric,
                        strict = strict,
                        min.bases = min.bases,
                        min.aligned.frac = min.aligned.frac,
                        isize.diff = isize.diff,
                        ...)

    if (all.reads) {
        window.reads.dt[, supporting := qname %in% rc$qname]
        return(window.reads.dt)
    }
    return(as.data.table(rc))
}


