\name{matchSamples}
\alias{matchSamples}
\title{Matches count data samples to meta data samples.}
\usage{
matchSamples(raw, filt.meta)
}
\description{
Makes a subset of the gene x sample matrix to match sample present in a subset of the meta data.}
\arguments{
\item{raw}{\code{data.frame}, gene x count matrix}
\item{filt.meta}{\code{data.frame}, meta data}
}
\value{
Returns a \code{data.frame} containing the gene x sample matrix for all samples in filt.meta.
}
\examples{
# filt.raw <- matchSamples(raw, filt.meta)
}
