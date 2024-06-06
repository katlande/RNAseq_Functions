\name{homer_to_matrix}
\alias{homer_to_matrix}
\title{Format a Homer raw count file for DESeq2 analysis}
\usage{
homer_to_matrix(raw, cutchar = "_", instance="1", pool=FALSE)
}
\description{
Takes the raw count output from Homer and reformats it for downstream analysis.}
\arguments{
\item{raw}{\code{data.frame}, Homer output file}
\item{cutchar}{\code{character}, removes all text from column names after the number of instances specified by \code{instance}}
\item{instance}{\code{numeric}, instances of \code{cutchar} after which to remove text. E.g., if \code{cutchar="_"} and \code{instance=2}, the name "sample_1_text" would become "sample_1".}
\item{pool}{\code{boolean}, if \code{TRUE}, reads from genes with identical names will be summed and used as a single entry. Prevents duplicate row.names if multiple transcripts are counted separately.}
\item{round}{\code{boolean}, whether or not to round fractional values to integers. TRUE is recommended for raw reads, but FALSE is recommended for any sort of normalized values (TPMs, FPKMs, etc)}
}
\value{
Returns a gene x sample \code{data.frame}.
}
\examples{
# raw <- homer_to_matrix(raw)
}
