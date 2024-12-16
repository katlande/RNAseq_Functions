\name{writeTSV}
\alias{writeTSV}
\title{save data.frame to TSV.}
\usage{
writeTSV(df, dir, rows=F)
}
\description{
Saves df as a .tsv in dir.}
\arguments{
\item{df}{An object to save as a .tsv}
\item{dir}{path to write output file to.}
\item{rows}{\code{boolean}, whether the output should save the row.names.}
}
\value{
Saves a file.
}
\examples{
# writeTSV(df, "/path/to/output", rows=F)
}
