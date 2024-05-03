\name{makeComp}
\alias{makeComp}
\title{Pull differential expression data from a DESeq2 object.}
\usage{
makeComp(Dnorm, var, up, down)
}
\description{
Pulls differential expression statistics from a comparison in a \code{DESeqDataSet} object.}
\arguments{
\item{Dnorm}{\code{DESeqDataSet}, a \code{DESeq_and_Filt()} output}
\item{var}{\code{character}, comparison variable name, must be a column name in the meta data}
\item{up}{\code{character}, level of \code{var} upregulated in the comparison, must match a level of the \code{var} column in the meta data}
\item{down}{\code{character}, level of \code{var} downregulated in the comparison, must match a level of the \code{var} column in the meta data}
}
\value{
Returns a \code{data.frame} containing differential expression statistics for a comparison.
}
\examples{
# DE_stats <- makeComp(Dnorm, "Age", "old", "young")
}
