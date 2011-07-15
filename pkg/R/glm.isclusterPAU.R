glmAndZIP.iscluster<-function(stfdf, idx, idxorder, minDateCluster, maxDateCluster, fractpop, modelCluster){
# Args:
# stfdf: data
# idTime, idSpace: indexes corresponding to the time and locations inside the cluster
# Returns:

# Fit null model
d0<-stfdf@data
switch(modelCluster,
poisson={m0<-glm(     Observed ~ 1,   offset=log(Expected), data=d0, family=poisson())}, 
zip=    {m0<-zeroinfl(Observed ~ 1|1, offset=log(Expected), data=d0)})

# difLaux is always >= 0. In the first iteration (i = 1) ncluster<-1 and difL<-difLaux for model i = 1.
sizeCluster<- -1
difL<- -1 
idTime<- which( (time(stfdf@time) >= minDateCluster) &  (time(stfdf@time) <= maxDateCluster))
#idSpace

if(length(idx) == 0) {
print('length(idx)=0')
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
}else{
return(c(sizeCluster, minDateCluster, maxDateCluster, difL))
}}
}
return(c(sizeCluster, minDateCluster, maxDateCluster, difL))
}




SetVbleCluster<-function(stfdf,idTime,idSpace){
# Args:
# stfdf: data
# idTime, idSpace: indexes corresponding to the time and locations inside the cluster
# Returns:
# vbleCluster, vbleCluster[i] = 1 if position i corresponds to space and time inside the cluster, 0 otherwise
vbleCluster<-rep(0,length(stfdf[['Observed']]))
idSubsetTime<-stfdf[,idTime,drop=FALSE]@data$ID
idSubsetSpace<-stfdf[idSpace,drop=FALSE]@data$ID
vbleCluster[(stfdf$ID %in% idSubsetTime) & (stfdf$ID %in% idSubsetSpace)]<-1
return(vbleCluster)
}
