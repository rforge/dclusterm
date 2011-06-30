glm.iscluster<-function(stfdf, idx, idxorder, minDateCluster, maxDateCluster, fractpop){
# Args:
# stdf: data
# idTime, idSpace: indexes corresponding to the time and locations inside the cluster
# Returns:

# Fit null model
m0<-glm(Observed ~ offset(log(Expected)), data=stfdf@data, family=poisson())

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
m0$data$CLUSTER<-SetVbleCluster(stdf,idTime,idSpace)
# cluster size must be smaller than fractpop of the total
if((sum(m0$data$CLUSTER*stfdf[['Expected']])-fractpop*sum(stfdf[['Expected']]))<0){
m1<-glm(update.formula(m0$formula, .~.+CLUSTER), data=m0$data, family=m0$family)
# 2 is the index for CLUSTER coefficient
# (-2*logLik(m0)+2*logLik(m1))/2
difLaux<-ifelse(coef(m1)[2]>0, (deviance(m0)-deviance(m1))/2, 0)
if(difLaux > difL){
sizeCluster<-i
difL<-difLaux
}else{
return(c(sizeCluster, minDateCluster, maxDateCluster, difL))
}}
}
return(c(sizeCluster, minDateCluster, maxDateCluster, difL))
}




SetVbleCluster<-function(stdf,idTime,idSpace){
# Args:
# stdf: data
# idTime, idSpace: indexes corresponding to the time and locations inside the cluster
# Returns:
# vbleCluster, vbleCluster[i] = 1 if position i corresponds to space and time inside the cluster, 0 otherwise
vbleCluster<-rep(0,length(stfdf[['Observed']]))
idSubsetTime<-stfdf[,idTime,drop=FALSE]@data$ID
idSubsetSpace<-stfdf[idSpace,drop=FALSE]@data$ID
vbleCluster[(stfdf$ID %in% idSubsetTime) & (stfdf$ID %in% idSubsetSpace)]<-1
return(vbleCluster)
}
