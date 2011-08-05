\name{CalcStatClusterGivenCenter}
\alias{CalcStatClusterGivenCenter}
\title{Calls the function to obtain the cluster with the maximum log-likelihood...}
\usage{CalcStatClusterGivenCenter(point, stfdf, rr, minDateCluster,
    maxDateCluster, fractpop, modelType, modelFormula,
    modelFamilyGlmGlmer, modelDistZeroinfl, modelLinkZeroinfl)}
\description{Calls the function to obtain the cluster with the maximum log-likelihood
ratio of all the clusters with the same center and start and end dates.}
\details{This function orders the regions according to the distance to a given center
and selects the regions with distance to the center less than sqrt(rr).
Then it calls glmAndZIP.iscluster() to obtain the cluster with the maximum
log-likelihood ratio of all the clusters with the same center and start and
end dates, and where the maximum fraction of the total population inside the
cluster is less than fractpop.}
\value{vector containing the coordinates of the center, the size, the
start and end dates, and the log-likelihood ratio of the cluster with the maximum log-likelihood ratio.}
\arguments{\item{point}{vector with the coordinates of the center of the cluster.}
\item{stfdf}{spatio-temporal class object containing the data.}
\item{rr}{square of the maximum radius of the cluster.}
\item{minDateCluster}{start date of the cluster.}
\item{maxDateCluster}{end date of the cluster.}
\item{fractpop}{maximum fraction of the total population inside the cluster.}
\item{modelType}{character specification of the type of probability model used to fit the data.
This can be "glm" for generalized linear models (glm {stats}),
"glmer" for generalized linear mixed model (glmer {lme4}), or
"zeroinfl" for zero-inflated models (zeroinfl {pscl}).}
\item{modelFormula}{character specification of the symbolic description of the model.}
\item{modelFamilyGlmGlmer}{family function to be used in the model if modelType is "glm" or "glmer".}
\item{modelDistZeroinfl}{character specification of count model family if modelType is "zeroinfl".}
\item{modelLinkZeroinfl}{character specification of link function in the binary zero-inflation model if modelType is "zeroinfl".}}
