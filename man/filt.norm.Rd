\name{filt.norm}
\alias{filt.norm}
\title{Normalize raw counts.}
\usage{
filt.norm(D, thresh=5, return.tab=T)
}
\description{
Filters and normalizes raw count data with \code{DESeq2::estimateSizeFactors()}}
\arguments{
\item{D}{A \code{DESeqDataSet} object}
\item{thresh}{\code{numeric}, minimum counts/gene for gene inclusion}
\item{return.tab}{\code{boolean}, if \code{TRUE}, returns a \code{data.frame} or normalized counts. If \code{FALSE}, returns a \code{DESeqDataSet} object}
}
\value{
Returns either a \code{data.frame} or \code{DESeqDataSet} containing normalized count data
}
\examples{
# filt.norm(D)
}
