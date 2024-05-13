\name{cluster_mat}
\alias{cluster_mat}
\title{Clustering function from pheatmap.}
\usage{
cluster_mat(mat, distance="correlation", method="ward.D2")
}
\description{Uses pheatmap clustering methods, but plots a heatmap as a grob that can be more easily manipulated.}
\arguments{
\item{mat}{\code{matrix}, numeric matrix of input data, clustering always occurs on ROWS}
\item{distance}{\code{character}, can use any clustering distance available in \code{pheatmap} clustering_distance_rows: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'}
\item{method}{\code{character}, hclust method for clustering}
}
\value{
Returns an \code{hclust} object.
}
\examples{
# cluster_mat(mat)
# cluster_mat(mat, "euclidean")
}

