% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{salmon_quant}
\alias{salmon_quant}
\title{GC content and lengths of 2000 intergenic regions}
\format{
A list of three elements:
\describe{
\item{abundance}{A numeric matrix containing abundance (TPM) for each gene
of each sample}
\item{counts}{A numeric matrix containing read count (fraction) for each
gene of each sample}
\item{length}{A numeric matrix containing length (bp) for each gene of
each sample}
}
}
\usage{
data(salmon_quant)
}
\description{
GC content and lengths of 2000 human intergenic regions calculated using the
\code{\link[=calc_region_gc]{calc_region_gc()}} function.
}
\keyword{datasets}
