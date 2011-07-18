\name{SetVbleCluster}
\alias{SetVbleCluster}
\title{Constructs a variable that indicates the locations and times that pertain...}
\usage{SetVbleCluster(stfdf, idTime, idSpace)}
\description{Constructs a variable that indicates the locations and times that pertain
to a cluster.}
\details{This function constructs a variable that indicates the locations and times
that pertain to a cluster. Each position of the variable is equal to 1 if
it corresponds to a location and time inside the cluster, and 0 otherwise.
This is one of the explanatory variables used in the glmAndZIP.iscluster
function to model the observed cases.}
\value{vector with 1's or 0's that indicates the locations and times that
pertain to a cluster.}
\arguments{\item{stfdf}{spatio-temporal class object containing the data.}
\item{idTime}{vector with the indexes of the stfdf object corresponding to
the time inside the cluster.}
\item{idSpace}{vector with the indexes of the stfdf object corresponding to
the locations inside the cluster.}}
