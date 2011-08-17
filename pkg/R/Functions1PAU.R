##' Detects clusters and computes their significance.
##' 
##' Searches all possible clusters with start and end dates within minDateUser
##' and maxDateUser, so that the maximum fraction of the total population inside
##' the cluster is less than fractpop, and the maximum distance to the center is
##' less than radius.
##' The search can be done for spatial or spatio-temporal clusters.
##' The significance of the clusters is obtained with a Monte Carlo procedure or
##' based on the chi-square distribution.
##'
##' @param stfdf spatio-temporal class object containing the data. See
##' STFDF-class {spacetime} for details. It contains an object of class
##' Spatial with the coordinates, a POSIXct object with the time, and a
##' data.frame with vectors Observed, Expected and potential covariates in each location and time.
##' @param thegrid two-columns matrix containing the points of the grid to be
##' used. If it is null, a rectangular grid is built.
##' @param radius maximum radius of the clusters.
##' @param step step of the thegrid built.
##' @param fractpop maximum fraction of the total population inside the cluster.
##' @param alpha significance level used to determine the existence of clusters.
##' @param typeCluster type of clusters to be detected. "ST" for spatio-temporal
##' or "S" spatial clusters.
##' @param minDateUser start date of the clusters.
##' @param maxDateUser end date of the clusters.
##' @param R If the cluster's significance is calculated based on the chi-square
##' distribution, R is NULL. If the cluster's significance is calculated using a
##' Monte Carlo procedure, R represents the number replicates under the null hypothesis.
##' @param numCPUS Number of cpus used when using snowfall to run the method.
##' If snowfall is not used numCPUS is NULL.
##' @param model0 Initial model (including covariates).
##' This can be "glm" for generalized linear models (glm {stats}),
##' "glmer" for generalized linear mixed model (glmer {lme4}), or
##' "zeroinfl" for zero-inflated models (zeroinfl {pscl}).
##'
##' @return data frame with information of the detected clusters ordered by its
##' log-likelihood ratio value. Each row represents the information of one of
##' the clusters. It contains the coordinates of the center, the size, the start
##' and end dates, the log-likelihood ratio, a boolean indicating if it is a
##' cluster (TRUE in all cases), and the p-value of the cluster.
##'
DetectClustersModel<-function(stfdf, thegrid=NULL, radius=Inf, step=NULL, fractpop, alpha,
typeCluster, minDateUser=min(time(stfdf@time)), maxDateUser=max(time(stfdf@time)), R=NULL, numCPUS=NULL, model0){

# Create column with ID. Unique identifier
stfdf[['ID']]<-1:nrow(stfdf@data)

sortDates<-sort(unique(time(stfdf@time)))

# Check minDateUser and maxDateUser make sense
minDateData<-min(time(stfdf@time))
maxDateData<-max(time(stfdf@time))
if(minDateUser > maxDateUser){
print('Error: cluster minimum date is greater than cluster maximum date')
return(0)
}
if(minDateUser > maxDateData){
print('Error: cluster minimum date is greater than data set maximum date')
return(0)
}
if(maxDateUser < minDateData){
print('Error: cluster maximum date is smaller than data set minimum date')
return(0)
}
# Closest dates to minDateUser and maxDateUser
idMinDateCluster<-min(which(sortDates >= minDateUser))
idMaxDateCluster<-max(which(sortDates <= maxDateUser))


# If grid is null, create a new grid
if(is.null(thegrid)){
CreateGridDClusterm(stfdf, radius, step)
}

# Radius
rr<-radius*radius

# Init snowfall
if(!is.null(numCPUS)){
sfInit( parallel=TRUE, cpus=numCPUS )
sfLibrary(spdep)
sfLibrary(splancs)
sfLibrary(spacetime)
sfLibrary(DCluster)
sfLibrary(pscl)
sfLibrary(DClusterm)
#sfSource("R/Functions1PAU.R")
#sfSource("R/Functions2PAU.R")
#sfSource("R/glm.isclusterPAU.R")
#sfSource("R/knutils.R")
}

# Statistic of each cluster
statsAllClusters<-CalcStatsAllClusters(thegrid, CalcStatClusterGivenCenter, stfdf, rr,
typeCluster, sortDates, idMinDateCluster, idMaxDateCluster, fractpop, model0, 
  numCPUS)

# Remove rows where sizeCluster == -1
idRemove<-which(statsAllClusters$sizeCluster == -1)
if(length(idRemove)>0){
statsAllClusters<-statsAllClusters[-idRemove, ]
}

# If there are no clusters return "No clusters found"
if(dim(statsAllClusters)[1] == 0){
print("No clusters found")
return("No clusters found")
}

# p-value of each cluster
vecpvalue<-matrix(NA,dim(statsAllClusters)[1],1)
veccluster<-matrix(NA,dim(statsAllClusters)[1],1)

##############################################################################################

# 1. p-value without Monte Carlo
if(is.null(R)){
for(i in 1:nrow(statsAllClusters)){
vecpvalue[i]<- 1-pchisq(2*statsAllClusters$statistic[i], 1)
veccluster[i]<-vecpvalue[i]<alpha
}}else{

# 2. p-value with Monte Carlo
maxStatisticRReplicas<-matrix(NA,R,1)
stfdfMC<-stfdf

for(i in 1:R){
# Generate data set under H_0
#stfdfMC$Observed<-rpois(length(stfdf$Observed), lambda = stfdf$Expected)

obslab<-as.character(formula(model0))[2]

stfdfMC[[obslab]]<-rpois(nrow(stfdf@data), lambda = fitted(model0))
# Statistic of each cluster
statsAllClustersMC<-CalcStatsAllClusters(thegrid, CalcStatClusterGivenCenter, stfdfMC, rr,
typeCluster, sortDates, idMinDateCluster, idMaxDateCluster, fractpop, model0,
  numCPUS)
maxStatisticRReplicas[i]<-max(statsAllClustersMC$statistic)
print(paste("replica",i))
}

# p-value according to rank
for(i in 1:(dim(statsAllClusters)[1])){
vecpvalue[i]<-(sum(maxStatisticRReplicas > statsAllClusters$statistic[i]) + 1)/(R + 1)
veccluster[i]<-vecpvalue[i]<alpha
}}

##############################################################################################

# End snowfall
if(!is.null(numCPUS)){
sfStop()
}

statsAllClusters<-cbind(statsAllClusters,veccluster,vecpvalue)
names(statsAllClusters)<-c("x", "y", "size", "minDateCluster", "maxDateCluster", "statistic", "cluster", "pvalue")

#print(statsAllClusters[rev(order(statsAllClusters$statistic)), ])
# Selection of significant clusters
statsAllClusters<-statsAllClusters[statsAllClusters$pvalue < alpha, ]

# If there are no clusters return "No clusters found"
if(dim(statsAllClusters)[1] == 0){
print(paste("No significant clusters found with alpha =", alpha))
return("No clusters found")
}

# Ordered results by statistic value
statsAllClusters<-statsAllClusters[rev(order(statsAllClusters$statistic)), ]

# Return
statsAllClusters$minDateCluster<-as.POSIXct(statsAllClusters$minDateCluster, origin="1970-01-01", tz="GMT")
statsAllClusters$maxDateCluster<-as.POSIXct(statsAllClusters$maxDateCluster, origin="1970-01-01", tz="GMT")
return(statsAllClusters)
}





##' Removes the overlapping clusters.
##' 
##' Function DetectClustersModel() detects duplicated clusters.
##' This function reduces the number of clusters by removing the overlapping
##' clusters.
##'
##' @param stfdf spatio-temporal class object containing the data.
##' @param statsAllClusters data frame with information of the detected
##' clusters obtained with DetectClustersModel().
##'
##' @return data frame with the same information than statsAllClusters but only
##' for clusters that do not overlap.
##'
SelectStatsAllClustersNoOverlap<-function(stfdf, statsAllClusters){
# statsAllClusters is ordered by statistic value
coordx<-as.data.frame(coordinates(stfdf@sp))[['x']]
coordy<-as.data.frame(coordinates(stfdf@sp))[['y']]
idSpaceAllClustersNoOverlap<-NULL
idTimeAllClustersNoOverlap<-NULL
statsAllClustersNoOverlap<-NULL

for(i in 1:(dim(statsAllClusters)[1])){
xd<-(coordx-statsAllClusters$x[i])
yd<-(coordy-statsAllClusters$y[i])
dist<-xd*xd+yd*yd
idSpaceOneCluster<-order(dist)[1:statsAllClusters$sizeCluster[i]]
idTimeOneCluster<-which(time(stfdf@time) >= statsAllClusters$minDateCluster[i] & time(stfdf@time) <= statsAllClusters$maxDateCluster[i])
if(sum(idSpaceOneCluster %in% idSpaceAllClustersNoOverlap) ==0 || sum(idTimeOneCluster %in% idTimeAllClustersNoOverlap) == 0){
statsAllClustersNoOverlap<-rbind(statsAllClustersNoOverlap, statsAllClusters[i, ])
idSpaceAllClustersNoOverlap<-c(idSpaceAllClustersNoOverlap, idSpaceOneCluster)
idTimeAllClustersNoOverlap<-c(idTimeAllClustersNoOverlap, idTimeOneCluster)
}}
return(statsAllClustersNoOverlap)
}








