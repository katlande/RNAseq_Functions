\name{display_venn}
\alias{display_venn}
\title{A pipe to VennDiagram.}
\usage{
display_venn(x, ...)
}
\description{
A pipe to the \code{VennDiagram::venn.diagram} function that draws the output using the \code{grid} package.}
\arguments{
\item{x}{A list of character vectors to overlap}
\item{...}{Additional parameters from \code{VennDiagram::venn.diagram}}
}
\value{
Returns a plotted venn diagram.
}
\examples{
# display_venn(list(genes_A, genes_B), 
#             category.names = c("Genes A" , "Genes B"), 
#             label=TRUE)
}
