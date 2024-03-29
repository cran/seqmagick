% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aln.R
\name{bs_aln}
\alias{bs_aln}
\title{bs_aln}
\usage{
bs_aln(x, method = "muscle", ...)
}
\arguments{
\item{x}{XStringSet object}

\item{method}{alignment method}

\item{...}{additional parameter}
}
\value{
aligned sequences, XStringSet object
}
\description{
sequence alignment
}
\examples{
\dontrun{
fa_file <- system.file("extdata/HA.fas", package="seqmagick")
x <- fa_read(fa_file)
bs_aln(x)
}
}
\author{
Guangchuang Yu
}
