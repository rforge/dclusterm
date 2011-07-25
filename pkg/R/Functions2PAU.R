##' Creates grid over the study area.
##' 
##' If the argument thegrid of DetectClustersModel() is null, this function is
##' used to create a rectangular grid with a given step.
##' If step is NULL the step used is equal to 0.2*radius.
##' The grid contains the coordinates of the centers of the clusters explored.
##'
##' @param stfdf spatio-temporal class object containing the data.
##' @param radius maximum radius of the clusters.
##' @param step step of the grid.
##'
##' @return two columns matrix where each row represents a point of the grid.
CreateGridDClusterm<-function(stfdf, radius, step){
# Return: thegrid
if(is.null(step)){
step<-.2*radius
}
coordx<-as.data.frame(coordinates(stfdf@sp))[['x']]
coordy<-as.data.frame(coordinates(stfdf@sp))[['y']]
xgrid<-seq(min(coordx), max(coordx), by=step)
ygrid<-seq(min(coordy), max(coordy), by=step)
xlen<-length(xgrid)
ylen<-length(ygrid)
npoints<-xlen*ylen
thegrid<-matrix(rep(NA, 2*npoints) , ncol=2)
thegrid[,1]<-rep(xgrid, times=ylen)
thegrid[,2]<-rep(ygrid, each=xlen)
return(thegrid)
}





##' Obtains clusters with the maximum log-likelihood ratio for each center and
##' start and end dates.
##' 
##' This function explores all possible clusters changing their center and start
##' and end dates. For each center and time periods, it obtains the cluster with
##' the maximum log-likelihood ratio so that the maximum fraction of the total
##' population inside the cluster is less than fractpop, and the maximum
##' distance to the center is less than radius.
##'
##' @param thegrid grid with the coordinates of the centers of the clusters
##' explored.
##' @param CalcStatClusterGivenCenter function to obtain the cluster with the
##' maximum log-likelihood ratio of all the clusters with the same center and
##' start and end dates
##' @param stfdf spatio-temporal class object containing the data.
##' @param rr square of the maximum radius of the cluster. 
##' @param typeCluster type of clusters to be detected. "ST" for spatio-temporal
##' clusters or "S" spatial clusters.
##' @param sortDates sorted vector of the times where disease cases occurred.
##' @param idMinDateCluster index of the closest date to the start date of the
##' cluster in the vector sortDates
##' @param idMaxDateCluster index of the closest date to the end date of the
##' cluster in the vector sortDates
##' @param fractpop maximum fraction of the total population inside the cluster.
##' @param modelCluster type of probability model used to fit the data. If
##' "poisson" generalized linear models with poisson family are used (glm {stats}).
##' If "zip" zero-inflated models are used (zeroinfl {pscl}).
##' @param numCPUS Number of cpus used when using snowfall to run the method.
##' If snowfall is not used numCPUS is NULL.
##'
##' @return data frame with information of the clusters with the maximum
##' log-likelihood ratio for each center and start and end dates. It contains the coordinates of the center, the size, the start and end dates, and the log-likelihood ratio of each of the clusters.
CalcStatsAllClusters<-function(thegrid, CalcStatClusterGivenCenter, stfdf, rr,
typeCluster, sortDates, idMinDateCluster, idMaxDateCluster, fractpop, modelCluster, numCPUS){
# Temporal dimension here, spatial dimension inside glmAndZIP.iscluster

if(typeCluster == "ST"){
statsAllClusters<-NULL
for (i in idMinDateCluster:idMaxDateCluster){
for (j in i: idMaxDateCluster){
if(is.null(numCPUS)){
statClusterGivenCenter<-apply(thegrid, 1, CalcStatClusterGivenCenter, stfdf, rr,
minDateCluster=sortDates[i], maxDateCluster=sortDates[j], fractpop, modelCluster)
}else{
statClusterGivenCenter<-sfApply(thegrid, 1, CalcStatClusterGivenCenter, stfdf, rr,
minDateCluster=sortDates[i], maxDateCluster=sortDates[j], fractpop, modelCluster)
}
statsAllClusters<-rbind(statsAllClusters,t(statClusterGivenCenter))
print(c(i,j))
}}}

if(typeCluster == "S"){
i<-idMinDateCluster
j<-idMaxDateCluster
if(is.null(numCPUS)){
statsAllClusters<-apply(thegrid, 1, CalcStatClusterGivenCenter, stfdf, rr,
minDateCluster=sortDates[i], maxDateCluster=sortDates[j], fractpop, modelCluster)
}else{
statsAllClusters<-sfApply(thegrid, 1, CalcStatClusterGivenCenter, stfdf, rr,
minDateCluster=sortDates[i], maxDateCluster=sortDates[j], fractpop, modelCluster)
}
statsAllClusters<-t(statsAllClusters)
print(c(i,j))
}
colnames(statsAllClusters)<-c("x", "y", "sizeCluster", "minDateCluster", "maxDateCluster", "statistic")
return(as.data.frame(statsAllClusters))
}





##' Calls the function to obtain the cluster with the maximum log-likelihood
##' ratio of all the clusters with the same center and start and end dates.
##' 
##' This function orders the regions according to the distance to a given center
##' and selects the regions with distance to the center less than sqrt(rr).
##' Then it calls glmAndZIP.iscluster() to obtain the cluster with the maximum
##' log-likelihood ratio of all the clusters with the same center and start and
##' end dates, and where the maximum fraction of the total population inside the
##' cluster is less than fractpop.
##'
##' @param point vector with the coordinates of the center of the cluster.
##' @param stfdf spatio-temporal class object containing the data.
##' @param rr square of the maximum radius of the cluster. 
##' @param minDateCluster start date of the cluster.
##' @param maxDateCluster end date of the cluster.
##' @param fractpop maximum fraction of the total population inside the cluster.
##' @param modelCluster type of probability model used to fit the data. If
##' "poisson" generalized linear models with poisson family are used (glm {stats}).
##' If "zip" zero-inflated models are used (zeroinfl {pscl}).
##'
##' @return vector containing the coordinates of the center, the size, the
##' start and end dates, and the log-likelihood ratio of the cluster with the maximum log-likelihood ratio.
CalcStatClusterGivenCenter<-function(point, stfdf, rr, minDateCluster, maxDateCluster, fractpop, modelCluster){
coordx<-as.data.frame(coordinates(stfdf@sp))[['x']]
coordy<-as.data.frame(coordinates(stfdf@sp))[['y']]
xd<-(coordx-point[1])
yd<-(coordy-point[2])
dist<-xd*xd+yd*yd
#
idx<-(dist <= rr)
idxorder<-order(dist)

# Only the regions with distance less than radius can be part of the cluster
idxorder<-idxorder[idx[idxorder]]

cl<-glmAndZIP.iscluster(stfdf=stfdf, idxorder=idxorder, minDateCluster, maxDateCluster, fractpop, modelCluster)
return(c(point, cl))
}


