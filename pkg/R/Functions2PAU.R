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


CalcStatsAllClusters<-function(thegrid, CalcStatClusterGivenCenter, stfdf, rr, iscluster,
typeCluster, sortDates, idMinDateCluster, idMaxDateCluster, fractpop){
# Return
# statsAllClusters
# Temporal dimension here, spatial dimension inside opgamModel

if(typeCluster == "ST"){
statsAllClusters<-NULL
for (i in idMinDateCluster:idMaxDateCluster){
for (j in i: idMaxDateCluster){
statClusterGivenCenter<-apply(thegrid, 1, CalcStatClusterGivenCenter, stfdf, rr, iscluster,
minDateCluster=sortDates[i], maxDateCluster=sortDates[j], fractpop)
statsAllClusters<-rbind(statsAllClusters,t(statClusterGivenCenter))
print(c(i,j))
}}}

if(typeCluster == "S"){
i<-idMinDateCluster
j<-idMaxDateCluster
statsAllClusters<-apply(thegrid, 1, CalcStatClusterGivenCenter, stfdf, rr, iscluster,
minDateCluster=sortDates[i], maxDateCluster=sortDates[j], fractpop)
statsAllClusters<-t(statsAllClusters)
print(c(i,j))
}
colnames(statsAllClusters)<-c("x", "y", "sizeCluster", "minDateCluster", "maxDateCluster", "statistic")
return(as.data.frame(statsAllClusters))
}


CalcStatClusterGivenCenter<-function(point, stfdf, rr, iscluster, minDateCluster, maxDateCluster, fractpop){
coordx<-as.data.frame(coordinates(stfdf@sp))[['x']]
coordy<-as.data.frame(coordinates(stfdf@sp))[['y']]
xd<-(coordx-point[1])
yd<-(coordy-point[2])
dist<-xd*xd+yd*yd
#
idx<-(dist <= rr)
idxorder<-order(dist)
cl<-iscluster(stfdf=stfdf, idx=idx, idxorder=idxorder, minDateCluster, maxDateCluster, fractpop)
return(c(point, cl))
}