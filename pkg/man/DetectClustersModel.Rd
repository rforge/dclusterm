\name{DetectClustersModel}
\alias{DetectClustersModel}
\title{Detects clusters and computes their significance.}
\usage{DetectClustersModel(stfdf, thegrid, radius=Inf, step, fractpop, alpha,
    typeCluster, minDateUser=min(time(stfdf@time)),
    maxDateUser=max(time(stfdf@time)), R, numCPUS, modelType="glm",
    modelFormula="Observed ~ 1", modelFamilyGlmGlmer=poisson(link =
    "log"), modelDistZeroinfl="poisson", modelLinkZeroinfl="logit")}
\description{Detects clusters and computes their significance.}
\details{Searches all possible clusters with start and end dates within minDateUser
and maxDateUser, so that the maximum fraction of the total population inside
the cluster is less than fractpop, and the maximum distance to the center is
less than radius.
The search can be done for spatial or spatio-temporal clusters.
The significance of the clusters is obtained with a Monte Carlo procedure or
based on the chi-square distribution.}
\value{data frame with information of the detected clusters ordered by its
log-likelihood ratio value. Each row represents the information of one of
the clusters. It contains the coordinates of the center, the size, the start
and end dates, the log-likelihood ratio, a boolean indicating if it is a
cluster (TRUE in all cases), and the p-value of the cluster.}
\arguments{\item{stfdf}{spatio-temporal class object containing the data. See
STFDF-class {spacetime} for details. It contains an object of class
Spatial with the coordinates, a POSIXct object with the time, and a
data.frame with vectors Observed, Expected and potential covariates in each location and time.}
\item{thegrid}{two-columns matrix containing the points of the grid to be
used. If it is null, a rectangular grid is built.}
\item{radius}{maximum radius of the clusters.}
\item{step}{step of the thegrid built.}
\item{fractpop}{maximum fraction of the total population inside the cluster.}
\item{alpha}{significance level used to determine the existence of clusters.}
\item{typeCluster}{type of clusters to be detected. "ST" for spatio-temporal
or "S" spatial clusters.}
\item{minDateUser}{start date of the clusters.}
\item{maxDateUser}{end date of the clusters.}
\item{R}{If the cluster's significance is calculated based on the chi-square
distribution, R is NULL. If the cluster's significance is calculated using a
Monte Carlo procedure, R represents the number replicates under the null hypothesis.}
\item{numCPUS}{Number of cpus used when using snowfall to run the method.
If snowfall is not used numCPUS is NULL.}
\item{modelType}{character specification of the type of probability model used to fit the data.
This can be "glm" for generalized linear models (glm {stats}),
"glmer" for generalized linear mixed model (glmer {lme4}), or
"zeroinfl" for zero-inflated models (zeroinfl {pscl}).}
\item{modelFormula}{character specification of the symbolic description of the model.}
\item{modelFamilyGlmGlmer}{family function to be used in the model if modelType is "glm" or "glmer".}
\item{modelDistZeroinfl}{character specification of count model family if modelType is "zeroinfl".}
\item{modelLinkZeroinfl}{character specification of link function in the binary zero-inflation model if modelType is "zeroinfl".}}
