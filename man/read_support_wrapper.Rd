% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contigsupport.R
\name{read_support_wrapper}
\alias{read_support_wrapper}
\title{read_support_wrapper}
\usage{
read_support_wrapper(
  le.dt = data.table(),
  reads.dt = data.table(),
  contigs.dt = data.table(),
  id = "",
  pad = 5000,
  ref.bwa = NULL,
  ref.fasta = NULL,
  bowtie = FALSE,
  outdir = ".",
  verbose = FALSE
)
}
\arguments{
\item{le.dt}{(full output from prep_loose_ends)}

\item{reads.dt}{(full output from loose_reads_wrapper)}

\item{contigs.dt}{($filtered.contigs from call_loose_end)}

\item{id}{(character) sample id}

\item{pad}{(numeric) window for getting supporting reads}

\item{ref.bwa}{(human output from grab_ref_obj)}

\item{ref.fasta}{(character) path to fasta if using bowtie}

\item{bowtie}{(logical) default FALSE}

\item{outdir}{(character) dirname for temporary files}

\item{verbose}{(logical)}
}
\value{
data.table with columns
\itemize{
\item leix
\item loose.end
\item contig.qname
\item tumor.count
\item normal.count
\item total.count
\item tumor.frac
}
}
\description{
read_support_wrapper
}
