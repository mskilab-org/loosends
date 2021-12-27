#' @name process_loose_ends
#' @title process_loose_ends
#'
#' @param id (character) sample ID
#' @param ranges (GRanges)
#' @param tbam (character) path to tumor bam
#' @param nbam (character) path to normal bam
#' @param ref_dir (character) path to reference directory
#' @param outdir (character) directory for storing outputs
#' @param read_pad (numeric) default 5000
#' @param assembly_pad (numeric) default 1000
#' @param verbose (logical) default TRUE
#'
#' @return list with two items
#' - call
#' - contigs
#' - support
process_loose_ends = function(id = "",
                              ranges = GRanges(),
                              tbam = "/dev/null",
                              nbam = "/dev/null",
                              ref_dir = "/dev/null",
                              outdir = "./",
                              read_pad = 5000,
                              assembly_pad = 1000,
                              verbose = FALSE) {

    human.ref = paste0(ref_dir, "/human_g1k_v37_decoy.fasta")
    reads.dt = loosereads_wrapper(ranges = ranges,
                               tbam = tbam,
                               nbam = nbam,
                               ref = human.ref,
                               outdir = outdir,
                               pad = read_pad,
                               verbose = verbose)

    le.dt = prep_loose_ends(li = ranges, id = id)
    ref.obj = grab_ref_obj(ref.dir = ref_dir)

    calls = call_loose_end_wrapper(id = id,
                                   le.dt = le.dt,
                                   reads.dt = reads.dt,
                                   ref_obj = ref.obj,
                                   pad = assembly_pad,
                                   mix.tn = TRUE,
                                   max.reads = 5000,
                                   verbose = verbose)

    support = read_support_wrapper(le.dt = le.dt,
                                   reads.dt = reads.dt,
                                   contigs.dt = calls$filtered.contigs,
                                   id = id,
                                   ref.bwa = ref.obj$human,
                                   verbose = verbose)


    return(list(calls = calls$calls,
                contigs = calls$filtered.contigs,
                support = support))
}
