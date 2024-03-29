% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loosereads.R
\name{.sample.spec2}
\alias{.sample.spec2}
\title{.sample.spec2}
\usage{
.sample.spec2(bam, chrsub = TRUE, verbose = FALSE)
}
\arguments{
\item{bam}{path to BAM file}

\item{chrsub}{(logical) substitute chr header? default TRUE}

\item{verbose}{optional, default=FALSE}
}
\description{
loads reads and mates for a single sample (tumor or normal)
assumes that BAM has already been filtered and avoids slow lapply
}
