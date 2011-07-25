\name{CalcStatsAllClusters}
\alias{CalcStatsAllClusters}
\title{Obtains clusters with the maximum log-likelihood ratio for each center and...}
\usage{CalcStatsAllClusters(thegrid, CalcStatClusterGivenCenter, stfdf, rr,
    typeCluster, sortDates, idMinDateCluster, idMaxDateCluster,
    fractpop, modelCluster, numCPUS)}
\description{Obtains clusters with the maximum log-likelihood ratio for each center and
start and end dates.}
\details{This function explores all possible clusters changing their center and start
and end dates. For each center and time periods, it obtains the cluster with
the maximum log-likelihood ratio so that the maximum fraction of the total
population inside the cluster is less than fractpop, and the maximum
distance to the center is less than radius.}
\value{data frame with information of the clusters with the maximum
log-likelihood ratio for each center and start and end dates. It contains the coordinates of the center, the size, the start and end dates, and the log-likelihood ratio of each of the clusters.}
\arguments{\item{thegrid}{grid with the coordinates of the centers of the clusters
explored.}
\item{CalcStatClusterGivenCenter}{function to obtain the cluster with the
maximum log-likelihood ratio of all the clusters with the same center and
start and end dates}
\item{stfdf}{spatio-temporal class object containing the data.}
\item{rr}{square of the maximum radius of the cluster.}
\item{typeCluster}{type of clusters to be detected. "ST" for spatio-temporal
clusters or "S" spatial clusters.}
\item{sortDates}{sorted vector of the times where disease cases occurred.}
\item{idMinDateCluster}{index of the closest date to the start date of the
cluster in the vector sortDates}
\item{idMaxDateCluster}{index of the closest date to the end date of the
cluster in the vector sortDates}
\item{fractpop}{maximum fraction of the total population inside the cluster.}
\item{modelCluster}{type of probability model used to fit the data. If
"poisson" generalized linear models with poisson family are used (glm {stats}).
If "zip" zero-inflated models are used (zeroinfl {pscl}).}
\item{numCPUS}{Number of cpus used when using snowfall to run the method.
If snowfall is not used numCPUS is NULL.}}
