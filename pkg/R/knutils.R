##' Make binary variables from the clusters.
##'
##' @param d whatever d is
##' @param knresults whatever knresults is
##'
##' @return whatever the return is
knbinary<-function(d, knresults){
clusters<-get.knclusters(d, knresults)
res<-lapply(clusters, function(X, n){
v<-rep(0,n)
v[X]<-1
return(v)}, n=dim(d)[1])

res<-data.frame(matrix(unlist(res), nrow=dim(d)[1]))
names(res)<-paste("CL", 1:length(clusters), sep="")
return(res)
}

##' Merge clusters so that they are identifed as factors.
##' Clusters should not overlap.
##'
##' @param d whatever d is
##' @param knresults whatever knresults is
##'
##' @return whatever the return is
mergeknclusters<-function(d, knresults){
n<-nrow(knresults)
knbin<-as.matrix(knbinary(d, knresults))	
res<-as.factor(knbin%*%matrix(1:n))
levels(res)<-c("NOCL", paste("CL", 1:n, sep="") )
return(res)
}

##' kn2SPDF.
##'
##' @param knresults whatever knresults is
##'
##' @return whatever the return is
kn2SPDF<-function(knresults){
clustercl<-SpatialPoints(knresults[,c("x", "y")] )
clustercl<-SpatialPointsDataFrame(clustercl, knresults)
}
