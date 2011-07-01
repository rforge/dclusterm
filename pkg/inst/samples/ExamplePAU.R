setwd("C:/GSOC11/GSOCJune21")

library(spdep)
library(RColorBrewer)
library(splancs)
library(spacetime)
library(DCluster)
library(pscl)


load("NY8.RData")
source("Functions1PAU.R")
source("Functions2PAU.R")
source("glm.isclusterPAU.R")
source("knutils.R")

# Calculate SMR
NY8$Exp<- NY8$POP8 * sum(NY8$Cases)/sum(NY8$POP8)
NY8$Observed<-NY8$Cases
NY8$Expected<-NY8$Exp
NY8$SMR<-NY8$Observed/NY8$Expected

# Add coordinates to the data.frame
xy<-coordinates(NY8)
NY8$x<-xy[,1]
NY8$y<-xy[,2]
names(xy)<-c("x", "y")

# Example 1
# STFDF object: sp, time and mydata
sp<-SpatialPoints(cbind(x = NY8$x[1:20], y = NY8$y[1:20]))
# time must be ordered
time<-as.POSIXct(strptime(c("1972-01-01", "1974-01-01", "1979-01-01", "1981-01-01",  "1983-01-01"), "%Y-%m-%d"), tz = "GMT")
mydata<-data.frame(
Observed = c(NY8$Cases[1:100]),
Expected = c(NY8$Exp[1:100]),
SMR      = c(NY8$SMR[1:100]))
stfdf = STFDF(sp, time, mydata)


# Example 2
# STFDF object: sp, time and mydata
sp<-SpatialPoints(cbind(x = NY8$x, y = NY8$y))
# time must be ordered
time<-as.POSIXct(strptime(c("1972-01-01"), "%Y-%m-%d"), tz = "GMT")
mydata<-data.frame(
Observed = c(NY8$Cases),
Expected = c(NY8$Exp),
SMR      = c(NY8$SMR))
stfdf = STFDF(sp, time, mydata)


# Create column with ID. Unique identifier
stfdf[['ID']]<-1:length(stfdf[['Observed']])

# Call method to detect clusters
# typeCluster="ST" (Spatio-temporal) or "S" (Spatial)
# modelCluster="poisson" (glm family poisson) or modelCluster="zip" (zeroinfl)

# if modelCluster="zip", stfdf$Observed<-round(stfdf$Observed)

statsAllClusters<-DetectClustersModel(stfdf=stfdf, thegrid=as.data.frame(stfdf@sp), radius=Inf, step=NULL, fractpop=0.15, alpha=0.05,
typeCluster="S", minDateUser=time(stfdf@time)[1], maxDateUser=time(stfdf@time)[1], modelCluster="zip")
statsAllClusters

# Select clusters that do not overlap
statsAllClustersNoOverlap<-SelectStatsAllClustersNoOverlap(stfdf,statsAllClusters)
statsAllClustersNoOverlap

# Plot detected clusters
knslim<-statsAllClustersNoOverlap[1, ]
knbinslim<-as.data.frame(knbinary(as(NY8, "data.frame"), knslim))
names(knbinslim)<-paste("CL", 1:ncol(knbinslim), sep="")

NY8test1<-NY8
NY8test1@data<-cbind(NY8test1@data, knbinslim)
NY8test1$clusters<-mergeknclusters(as(NY8test1, "data.frame"), knslim)

pal<-c("white", brewer.pal(12, "Set3"))
spplot(NY8test1, "clusters",  col="gray",  col.regions=pal[1:nlevels(NY8test1$clusters)],
sp.layout=list(list("sp.points", kn2SPDF(knslim), pch=19), list("sp.text", as.matrix(knslim[,1:2]), levels(NY8test1$clusters)[-1],
font=2, cex=1, col="black")))



