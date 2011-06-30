

#This function makes binary variables from the clusters
knbinary<-function(d, knresults)
{
	clusters<-get.knclusters(d, knresults)

	res<-lapply(clusters, function(X, n)
		{
			v<-rep(0,n)
			v[X]<-1

			return(v)
		}

		, n=dim(d)[1])

	res<-data.frame(matrix(unlist(res), nrow=dim(d)[1]))
	names(res)<-paste("CL", 1:length(clusters), sep="")

	return(res)
}


#This function slims the number of clusters down.
#The spatial scan statistic is known to detect duplicated
#clusters. This function aims to reduce the number of clusters
#by removing duplicated and overlapping clusters.
#
#The main criteria to choose the 'primary' clusters is the
#is the likelihood ratio


slimknclusters<-function(d, knresults, minsize=1)
{
	#Filter by minsize
	knresults<-knresults[which(knresults$size>=minsize),]

	#Ordering according to the test statistic
	idxcl<-rev(order(knresults$statistic))

	knbin<-knbinary(d, knresults)

	clusters<-c()
	while(length(idxcl)>0)
	{

		print(knresults[idxcl[1], ])

		cl<-idxcl[1]


		if(length(idxcl)>0)
		{
		res<-apply(as.matrix(knbin[,idxcl]), 2, 
			function(X, clbin){sum(X*clbin)}, 
			clbin=knbin[,cl])

	idxrem<-which(res>0)#Here is where we decide what clusters to remove
		idxcl<-idxcl[-c(1, idxrem)]
		}
		else
		{
		idxcl<-c()#idxcl[-c(1)]
		}

		clusters<-c(clusters, cl)
	}

	return(knresults[clusters,])
}


#Merge clusters so that they are identifed as factors
#Clusters should not overlap
mergeknclusters<-function(d, knresults)
{
	n<-nrow(knresults)
	knbin<-as.matrix(knbinary(d, knresults))
	
	res<-as.factor(knbin%*%matrix(1:n))
	levels(res)<-c("NOCL", paste("CL", 1:n, sep="") )

	return(res)


}


kn2SPDF<-function(knresults)
{
	clustercl<-SpatialPoints(knresults[,c("x", "y")] )
	clustercl<-SpatialPointsDataFrame(clustercl, knresults)
}


#
#Returns an index indicating how many 
#regions in clusters from  knresults1 overlap with cluster
#regions in knresults0
overlapindex<-function(d, knresults0, knresults1)
{
	knbin0<-knbinary(d, knresults0)
	knbin1<-knbinary(d, knresults1)
	
	res<-as.data.frame(matrix(NA, nrow=nrow(d), ncol=nrow(knresults0)))

	res<-apply(knbin0, 2, function(X){
			ov<-rep(0, length(X))

			ov[X==1]<-apply(knbin1[X==1,], 1,sum) 

			return(ov)
		})

	res<-as.data.frame(res)
		
	names(res)<-paste("IDXCL", 1:ncol(res), sep="")

	return(res)
}

#RETURN names of areas in cluster as a list
clusternames<-function(d, knresults, namescol)
{
	knbin<-knbinary(xx@data, knresults)
	apply(knbin, 2, function(X){
		d[which(X==1), namescol]
		})	
}

#Relevel clusters from mergeknclusters
relev<-function(x)
{
x<-as.factor(as.numeric(x)>1)
levels(x)<-c("", "CLUSTER")
x
}

