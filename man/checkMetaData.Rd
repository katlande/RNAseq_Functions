\name{checkMetaData}
\alias{checkMetaData}
\title{Checks if sample names in meta data and count data match.}
\usage{
checkMetaData(meta, counts)
}
\description{
Checks if meta data and count data match; reorders meta data when necessary, and warns users if files do not match.}
\arguments{
\item{meta}{A \code{data.frame} of meta data with sample names as \code{row.names}}
\item{counts}{A sample x gene \code{data.frame} with sample names as \code{col.names}}
}
\value{
Returns meta data ordered to match count data sample order.
}
\examples{
# meta <- checkMetaData(meta, counts)
}
