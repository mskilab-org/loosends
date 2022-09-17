#' @name find_telomeres2
#' @title find_telomeres2
#'
#' @param seq (character) vector of sequences
#' @param verbose (logical) default FALSE
#'
#' @return data.table with rows equal to the number of strings in seq and columns
#' - grtr_canonical
#' - grtr_noncanonical
#' - crtr_canonical
#' - crtr_noncanonical
#'
#' which are logical vectors showing whether there is a canonical/noncanonical GRTR/CRTR in each string in seq
find_telomeres2 = function(seq = character(), verbose = FALSE)
{
    if (!length(seq)) {
        return(data.table(grtr_canonical = logical(),
                          crtr_canonical = logical(),
                          grtr_noncanonical = logical(),
                          crtr_noncanonical = logical()))
    }

    if (verbose) {message("loading pdicts")}
    grtr.canonical.pdict = readRDS(system.file("extdata", "g.canonical.pdict.rds", package = "loosends"))
    crtr.canonical.pdict = readRDS(system.file("extdata", "c.canonical.pdict.rds", package = "loosends"))
    grtr.noncanonical.pdict = readRDS(system.file("extdata", "g.noncanonical.pdict.rds",
                                                  package = "loosends"))
    crtr.noncanonical.pdict = readRDS(system.file("extdata", "c.noncanonical.pdict.rds",
                                                  package = "loosends"))

    dstring = Biostrings::DNAStringSet(x = seq)

    ## check for perfect matches to telomeres
    if (verbose) { message("checking for exact telomeric matches") }
    grtr.canonical.res = Biostrings::vwhichPDict(pdict = grtr.canonical.pdict,
                                                 subject = dstring,
                                                 min.mismatch = 0,
                                                 max.mismatch = 0)
    crtr.canonical.res = Biostrings::vwhichPDict(pdict = crtr.canonical.pdict,
                                                 subject = dstring,
                                                 min.mismatch = 0,
                                                 max.mismatch = 0)
    grtr.noncanonical.res = Biostrings::vwhichPDict(pdict = grtr.noncanonical.pdict,
                                                    subject = dstring,
                                                    min.mismatch = 0,
                                                    max.mismatch = 0)
    crtr.noncanonical.res = Biostrings::vwhichPDict(pdict = crtr.noncanonical.pdict,
                                                    subject = dstring,
                                                    min.mismatch = 0,
                                                    max.mismatch = 0)

    return(data.table(grtr_canonical = base::lengths(grtr.canonical.res) > 0,
                      crtr_canonical = base::lengths(crtr.canonical.res) > 0,
                      grtr_noncanonical = base::lengths(grtr.noncanonical.res) > 0,
                      crtr_noncanonical = base::lengths(crtr.noncanonical.res) > 0))
                      
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

#' @name eighteenmer_from_fasta
#' @title eighteenmer_from_fasta
#'
#' @description
#'
#' Given a .fasta file of (forward strand) telomeric monomers, generate all possible eighteenmers
#'
#' @param fasta (character) path to ref fasta
#' @param rc (logical) reverse complement
#'
#' @return character vector of all possible eighteenmers
eighteenmer_from_fasta = function(fasta = system.file('extdata', 'telomeres.fa', package='loosends'),
                                  rc = FALSE)
{
    monomers = readDNAStringSet(filepath = fasta, format = "fasta")

    ## check validity of monomer lengths
    mlengths = base::lengths(monomers)
    if (!all(mlengths == 6)) {
        stop("not all monomers have length 6")
    }

    ## reverse complement monomers if applicable
    if (rc) {
        monomers = Biostrings::reverseComplement(x = monomers)
    }

    mstrings = unname(as.character(monomers))
    grid.dt = as.data.table(expand.grid(m1 = mstrings,
                                        m2 = mstrings,
                                        m3 = mstrings,
                                        m4 = mstrings,
                                        start.ix = 1:6))

    grid.dt[, seq := substr(paste0(m1, m2, m3, m4), start = start.ix, stop = start.ix + 17)]
    return(grid.dt[, seq])
}

#' @name eighteenmer
#' @title eighteenmer
#'
#' @description
#' 
#' generate a PDict of telomeric 18mers
#'
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

#' @name munch
#' @title munch
#'
#' @description
#' 
#' returns logical vector of length reads
#' indicating does read contain any match to query
#' searches for matches on same strand as input seq
#' 
#' @param reads GRanges or data.table containing $seq field with character sequence
#' @param query optional, PDict of motifs, default= 18mers of telomeric motifs (both strands)
munch = function(reads, query=NULL) {
    if(is.null(query)) query = eighteenmer()
    seq = reads$seq
    seq_frame_ref = DNAStringSet(seq)
    return(match.seq(query, seq_frame_ref))
}
