\name{add_col}
\alias{add_col}
\title{Adds a new column to a data.frame}
\usage{
CheckDepth(df, colname, colval)
}
\description{
Add a dataframe column.}
\arguments{
\item{df}{input data.frame}
\item{colname}{character, colname of new column.}
\item{colval}{a single value or a vector the length of nrow(df) containing the contents of the new column}
}
\value{
Returns a data.frame with ncol()+1 columns.
}
\examples{
# CheckDepth(df, "Set", "batch 1")
}
