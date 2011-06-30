DetectClustersModel<-function(stfdf, thegrid=NULL, radius=Inf, step=NULL,
fractpop, alpha, iscluster, typeCluster, minDateUser=NULL, maxDateUser=NULL){

idMinDateCluster<-NULL
idMaxDateCluster<-NULL
sortDates<-sort(unique(time(stfdf@time)))

if(typeCluster == "S"){
idMinDateCluster<-which.min(sortDates)
idMaxDateCluster<-which.max(sortDates)
}else{
if(typeCluster == "ST"){
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
}}

# If grid is null, create a new grid
if(is.null(thegrid)){
CreateGridDClusterm(stfdf, radius, step)
}

# Radius
rr<-radius*radius

# Statistic of each cluster
statsAllClusters<-CalcStatsAllClusters(thegrid, CalcStatClusterGivenCenter, stfdf, rr, iscluster,
typeCluster, sortDates, idMinDateCluster, idMaxDateCluster, fractpop)

# Remove rows where sizeCluster == -1
idRemove<-which(statsAllClusters$sizeCluster == -1)
if(length(idRemove)>0){
statsAllClusters<-statsAllClusters[-idRemove, ]
print("remove")
}

# If there are no clusters return "No clusters found"
if(dim(statsAllClusters)[1] == 0){
print("No clusters found")
return("No clusters found")
}

# p-value of each cluster
vecpvalue<-matrix(NA,dim(statsAllClusters)[1],1)
veccluster<-matrix(NA,dim(statsAllClusters)[1],1)
for(i in 1:(dim(statsAllClusters)[1])){
vecpvalue[i]<- 1-pchisq(2*statsAllClusters$statistic[i], 1)
veccluster[i]<-vecpvalue[i]<alpha
}
statsAllClusters<-cbind(statsAllClusters,veccluster,vecpvalue)
colnames(statsAllClusters)<-c("x", "y", "sizeCluster", "minDateCluster", "maxDateCluster", "statistic", "cluster", "pvalue")

# Selection of significant clusters
statsAllClusters<-statsAllClusters[statsAllClusters$pvalue < alpha, ]

# Ordered results by statistic value
statsAllClusters<-statsAllClusters[rev(order(statsAllClusters$statistic)), ]

# Return
statsAllClusters$minDateCluster<-as.POSIXct(statsAllClusters$minDateCluster, origin="1970-01-01", tz="GMT")
statsAllClusters$maxDateCluster<-as.POSIXct(statsAllClusters$maxDateCluster, origin="1970-01-01", tz="GMT")
return(statsAllClusters)
}




SelectStatsAllClustersNoOverlap<-function(stfdf,statsAllClusters){
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








