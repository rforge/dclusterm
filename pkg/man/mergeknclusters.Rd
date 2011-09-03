\name{mergeknclusters}
\alias{mergeknclusters}
\title{Merges clusters so that they are identifed as levels of a factor.}
\usage{mergeknclusters(datamap, knresults, indClustersPlot)}
\description{Merges clusters so that they are identifed as levels of a factor.}
\details{Given a data frame with clusters that do not overlap 
this function merges the clusters and construct a factor.
The levels of the factor are "NCL" if the polygon of the map is not
in any cluster, and "CLi" if the polygon i is in cluster i.}
\value{factor with levels that represent the clusters.}
\arguments{\item{datamap}{data of the SpatialPolygonsDataFrame with the polygons
of the map.}
\item{knresults}{Data frame with information of the detected clusters.
Each row represents the information of one of the clusters.
It contains the coordinates of the center, the size, the start
and end dates, the log-likelihood ratio, a boolean indicating if it is a
cluster (TRUE in all cases), and the p-value of the cluster.}
\item{indClustersPlot}{rows of knresults that denote the clusters to be plotted.}}
