##' Obtains the cluster with the maximum log-likelihood ratio of all the
##' clusters with the same center and start and end dates.
##' 
##' This function constructs all the clusters with start date equal to
##' minDateCluster, end date equal to maxDateCluster, and with center specified
##' by the first element of idxorder, so that the maximum fraction of the total
##' population inside the cluster is less than fractpop, and the maximum
##' distance to the center is less than radius.
##' For each one of these clusters, the log-likelihood ratio test statistic
##' for comparing the alternative model with the cluster versus the null model
##' of no clusters is calculated.
##' The cluster with maximum value of the log-likelihood ratio is returned.
##'
##' @param stfdf a spatio-temporal class object containing the data.
##' @param idxorder a permutation of the regions according to their distance to
##' the current center.
##' @param minDateCluster start date of the cluster.
##' @param maxDateCluster end date of the cluster.
##' @param fractpop maximum fraction of the total population inside the cluster.
##' @param modelCluster type of probability model used to fit the data. If
##' "poisson" generalized linear models with poisson family are used (glm {stats}). If "zip" zero-inflated models are used (zeroinfl {pscl}).
##'
##' @return vector containing the size, the start and end dates, and the
##' log-likelihood ratio of the cluster with the maximum log-likelihood ratio.
glmAndZIP.iscluster<-function(stfdf, idxorder, minDateCluster, maxDateCluster, fractpop, modelCluster){
# Fit null model
d0<-stfdf@data
switch(modelCluster,
poisson={m0<-glm(     Observed ~ 1,   offset=log(Expected), data=d0, family=poisson())}, 
zip=    {m0<-zeroinfl(Observed ~ 1|1, offset=log(Expected), data=d0)})

# difLaux is always >= 0. In the first iteration (i = 1) ncluster<-1 and difL<-difLaux for model i = 1.
sizeCluster<- -1
difL<- -1 
idTime<- which( (time(stfdf@time) >= minDateCluster) &  (time(stfdf@time) <= maxDateCluster))
# idTime, idSpace: indexes corresponding to the time and locations inside the cluster

if(length(idxorder) == 0) {
print('length(idxorder)=0')
return(c(sizeCluster, minDateCluster, maxDateCluster, difL))
}

for(i in 1:length(idxorder)){
idSpace<-idxorder[1:i]
d0$CLUSTER<-SetVbleCluster(stfdf,idTime,idSpace)

# cluster size must be smaller than fractpop of the total
# 2 is the index for CLUSTER coefficient

if((sum(d0$CLUSTER*stfdf[['Expected']])-fractpop*sum(stfdf[['Expected']]))<0){

switch(modelCluster,
poisson={m1<-glm(     Observed ~ CLUSTER,   offset=log(Expected), data=d0, family=m0$family)
difLaux<-ifelse(coef(m1)[2]>0, (deviance(m0)-deviance(m1))/2,0) }, 
zip    ={m1<-zeroinfl(Observed ~ CLUSTER|1, offset=log(Expected), data=d0)
difLaux<-ifelse(coef(m1)[2]>0, (-2*logLik(m0)+2*logLik(m1))/2, 0) })

if(difLaux > difL){
sizeCluster<-i
difL<-difLaux
}}else{
return(c(sizeCluster, minDateCluster, maxDateCluster, difL))
}
}
return(c(sizeCluster, minDateCluster, maxDateCluster, difL))
}





##' Constructs a variable that indicates the locations and times that pertain
##' to a cluster.
##' 
##' This function constructs a variable that indicates the locations and times
##' that pertain to a cluster. Each position of the variable is equal to 1 if
##' it corresponds to a location and time inside the cluster, and 0 otherwise.
##' This is one of the explanatory variables used in the glmAndZIP.iscluster
##' function to model the observed cases.
##'
##' @param stfdf spatio-temporal class object containing the data.
##' @param idTime vector with the indexes of the stfdf object corresponding to
##' the time inside the cluster.
##' @param idSpace vector with the indexes of the stfdf object corresponding to
##' the locations inside the cluster.
##'
##' @return vector with 1's or 0's that indicates the locations and times that
##' pertain to a cluster.
SetVbleCluster<-function(stfdf,idTime,idSpace){
# vbleCluster, vbleCluster[i] = 1 if position i corresponds to space and time inside the cluster, 0 otherwise
vbleCluster<-rep(0,length(stfdf[['Observed']]))
idSubsetTime<-stfdf[,idTime,drop=FALSE]@data$ID
idSubsetSpace<-stfdf[idSpace,drop=FALSE]@data$ID
vbleCluster[(stfdf$ID %in% idSubsetTime) & (stfdf$ID %in% idSubsetSpace)]<-1
return(vbleCluster)
}
