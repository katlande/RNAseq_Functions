\name{Heatmap_Grob}
\alias{Heatmap_Grob}
\title{Make a pheatmap as a grob.}
\usage{
Heatmap_Grob(mat, font_row=8, font_col=8, tilt=T, show_row=T, show_col=T, border="black",
clustDist_col="correlation", clustDist_row="correlation", clustMethod="ward.D2",
fill_colours=NULL, title=NULL, legend_title=NULL, scale_rows=T)
}
\description{Uses pheatmap clustering methods, but plots a heatmap as a grob that can be more easily manipulated.}
\arguments{
\item{mat}{\code{matrix}, numeric matrix of input data}
\item{font_row}{\code{numeric}, row label font size}
\item{font_col}{\code{numeric}, column label font size}
\item{tilt}{\code{boolean}, whether to tilt column labels}
\item{show_row}{\code{boolean}, whether to plot row labels}
\item{show_col}{\code{boolean}, whether to plot column labels}
\item{border}{\code{character}, colour for cell borders}
\item{clustDist_col}{\code{character}, pipe to \code{pheatmap} clustering_distance_cols}
\item{clustDist_row}{\code{character}, pipe to \code{pheatmap} clustering_distance_rows}
\item{clustMethod}{\code{character}, hclust method for clustering}
\item{fill_colours}{\code{character}, custom vector of colours to use for cell colours}
\item{title}{\code{character}, plot title}
\item{legend_title}{\code{character}, legend title}
\item{scale_rows}{\code{boolean}, whether to plot row-scaled Z-scores of matrix values}
}
\value{
Returns a \code{grob} containing a pheatmap object.
}
\examples{
# Heatmap_Grob(mat)
# Heatmap_Grob(mat, title="my heatmap", scale_rows=F, clustDist_col="euclidean", show_row=F)
}

