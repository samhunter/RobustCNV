% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/functions.R
\name{summarizeToGeneLevel}
\alias{summarizeToGeneLevel}
\title{Given a set of interval calls with gene information, use a set of rules to summarize them to
gene level calls.}
\usage{
summarizeToGeneLevel(interval.calls)
}
\arguments{
\item{intervals}{a dataframe with a definition of the intervals used in the capture experiment}

\item{gsegments}{dataframe returned from the callSegments() function.}

\item{normalized_tumor_coverage}{a vector returned from robustNorm() function.}
}
\value{
a dataframe combining intervals and calls
}
\description{
Given a set of interval calls with gene information, use a set of rules to summarize them to
gene level calls.
}

