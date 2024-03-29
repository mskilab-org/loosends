% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contig.support.R
\name{check_split_contig_support}
\alias{check_split_contig_support}
\title{check_split_contig_support}
\usage{
check_split_contig_support(calns, reads.dt, ref, seed.pad = 0, verbose = FALSE)
}
\arguments{
\item{calns}{(data.table) represents contig alignment}

\item{reads.dt}{(data.table) reads corresponding to the loose end associated with this contig}

\item{ref}{(BWA) reference to compare alignment against}

\item{seed.pad}{(numeric) number of base pairs around peak to search for supporting reads, default 0}

\item{verbose}{(logical) default FALSE}
}
\description{
Gets supporting reads for split contigs (aligning to two separate places in the reference)
Note that this is not appropriate for contigs that represent purely the mate sequence

Takes as input a contig alignment, reads, and a reference BWA object
}
