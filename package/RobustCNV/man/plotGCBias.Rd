% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plotGCBias}
\alias{plotGCBias}
\title{Plot GC bias plot}
\usage{
plotGCBias(sampleName, pGC, coverage)
}
\arguments{
\item{sampleName}{sample name}

\item{pGC}{vector of percent GC for each interval}

\item{coverage}{vector of per-interval coverage}
}
\value{
A gc-score value which is the sum of abs(predicted values).
}
\description{
Plot GC bias plot
}
\examples{
\dontrun{
plotGCBias(sampleName, pGC, coverage)
}
}

