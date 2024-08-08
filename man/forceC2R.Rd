\name{forceC2R}
\alias{forceC2R}
\title{Force tibble's column_to_rownames() when row.names are already present.}
\usage{
forceC2R(df, col, keep=F, colname="KEEP")
}
\description{
Run tibble's column_to_rownames on a dataframe with existing row.names.}
\arguments{
\item{df}{A \code{data.frame}}
\item{col}{\code{character}, the colname to set as row.names}
\item{keep}{\code{boolean}, whether or not to save the row.names in a column}
\item{colname}{\code{character}, if keep==T, the name of the new column}
}
\value{
Returns a \code{data.frame} with new row.names
}
\examples{
# forceC2R(df, "Treatment")
# forceC2R(df, "Treatment", T, "NewColumn")
}
