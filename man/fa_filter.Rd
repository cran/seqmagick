% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{fa_filter}
\alias{fa_filter}
\title{fa_filter}
\usage{
fa_filter(
  fasfile,
  pattern,
  by = "description",
  ignore.case = FALSE,
  outfile = NULL,
  type = "interleaved"
)
}
\arguments{
\item{fasfile}{input fasta file}

\item{pattern}{keyword for filter}

\item{by}{one of 'description' and 'sequence'}

\item{ignore.case}{logical}

\item{outfile}{output file}

\item{type}{one of 'interleaved' and 'sequential'}
}
\value{
BStringSet object
}
\description{
fasta filter by searching pattern
}
\author{
Guangchuang Yu
}
