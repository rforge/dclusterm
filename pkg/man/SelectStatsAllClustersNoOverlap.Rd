\name{SelectStatsAllClustersNoOverlap}
\alias{SelectStatsAllClustersNoOverlap}
\title{Removes the overlapping clusters.}
\usage{SelectStatsAllClustersNoOverlap(stfdf, statsAllClusters)}
\description{Removes the overlapping clusters.}
\details{Function DetectClustersModel() detects duplicated clusters.
This function reduces the number of clusters by removing the overlapping
clusters.}
\value{data frame with the same information than statsAllClusters but only
for clusters that do not overlap.}
\arguments{\item{stfdf}{spatio-temporal class object containing the data.}
\item{statsAllClusters}{data frame with information of the detected
clusters obtained with DetectClustersModel().}}
