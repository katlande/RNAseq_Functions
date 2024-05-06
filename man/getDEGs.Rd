\name{getDEGs}
\alias{getDEGs}
\title{Get DEGs.}
\usage{
getDEGs(results, mode="count", alpha=0.05, lfc=0.5, direction="all")
}
\description{
Get the number, names, or statistics of all DEGs from a \code{makeComp()} output \code{data.frame}.}
\arguments{
\item{results}{A \code{makeComp()} output \code{data.frame}}
\item{mode}{\code{character}, "count", "vector", or "table". Returns either the number of DEGs, a vector of all DEGs, or a filtered version of \code{results} with only DEGs}
\item{alpha}{\code{numeric}, pAdj threshold}
\item{lfc}{\code{numeric}, absolute log2foldchange threshold}
\item{direction}{\code{character}, "all", "up", or "down", whether to include all DEGs, upregulated DEGs only, or downregulated DEGs only. Only relevant for \code{mode == "vector"} or \code{mode == "table"}}
}
\value{
Returns a console \code{string}, a \code{vector}, or a \code{data.frame}, depending on \code{mode}.
}
\examples{
# getDEGs(results)
# getDEGs(results, mode="table", lfc=1, direction="up")
}
