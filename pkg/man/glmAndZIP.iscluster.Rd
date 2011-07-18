\name{glmAndZIP.iscluster}
\alias{glmAndZIP.iscluster}
\title{Obtains the cluster with the maximum log-likelihood ratio of all the...}
\usage{glmAndZIP.iscluster(stfdf, idxorder, minDateCluster, maxDateCluster,
    fractpop, modelCluster)}
\description{Obtains the cluster with the maximum log-likelihood ratio of all the
clusters with the same center and start and end dates.}
\details{This function constructs all the clusters with start date equal to
minDateCluster, end date equal to maxDateCluster, and with center specified
by the first element of idxorder, so that the maximum fraction of the total
population inside the cluster is less than fractpop, and the maximum
distance to the center is less than radius.
For each one of these clusters, the log-likelihood ratio test statistic
for comparing the alternative model with the cluster versus the null model
of no clusters is calculated.
The cluster with maximum value of the log-likelihood ratio is returned.}
\value{vector containing the size, the start and end dates, and the
log-likelihood ratio of the cluster with the maximum log-likelihood ratio.}
\arguments{\item{stfdf}{a spatio-temporal class object containing the data.}
\item{idxorder}{a permutation of the regions according to their distance to
the current center.}
\item{minDateCluster}{start date of the cluster.}
\item{maxDateCluster}{end date of the cluster.}
\item{fractpop}{maximum fraction of the total population inside the cluster.}
\item{modelCluster}{type of probability model used to fit the data. If
"poisson" generalized linear models with poisson family are used (glm {stats}). If "zip" zero-inflated models are used (zeroinfl {pscl}).}}
