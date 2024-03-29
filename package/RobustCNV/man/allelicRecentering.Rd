% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/functions.R
\name{allelicRecentering}
\alias{allelicRecentering}
\title{Some samples which are highly aneuploid are not well centered after normalization. This function
attempts to center them using allele fraction information to identify putatively diploid segments.}
\usage{
allelicRecentering(normCoverage, intervals, sampleAlleles)
}
\arguments{
\item{normCoverage}{normalized log2 ratios from robustNorm()}

\item{intervals}{dataframe defining the intervals}

\item{sampleAlleles}{dataframe with allele fractions and index}
}
\value{
a vector which has been re-centered.
}
\description{
Some samples which are highly aneuploid are not well centered after normalization. This function
attempts to center them using allele fraction information to identify putatively diploid segments.
}

