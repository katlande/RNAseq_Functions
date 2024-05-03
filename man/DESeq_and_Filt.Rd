\name{DESeq_and_Filt}
\alias{DESeq_and_Filt}
\title{Pipe to DESeq.}
\usage{
DESeq_and_Filt(D, filter_thresh, minSamples=3)
}
\description{
Filters data and runs DESeq2.}
\arguments{
\item{D}{A \code{DESeqDataSet} object}
\item{filter_thresh}{\code{numeric}, minimum total raw reads across all samples for inclusion}
\item{minSamples}{\code{numeric}, minimum number of samples with > 0 reads required for gene inclusion}
}
\value{
Returns a \code{DESeqDataSet} object.
}
\examples{
# DESeq_and_Filt(D, 5)
}
