library(cluster)
library(clusterSim)
library(vegan)
library(parallel)

data <- read.table("cor_mat_SparCC.out",header = T,row.names = 1)
pval <- read.table("pval_two_sided.txt",header=T,row.names = 1)

clus<-makeCluster(4)
clusterEvalQ(clus,library(vegan))
data.dist<-as.dist(1-data)
data.cluster<-hclust(data.dist,method="ward.D2")
# get.Fratio<-function(dist,cluster,maxgroupnum){
	# Fratio<-NULL
	# cluster.cut<-cutree(cluster,c(1:maxgroupnum))
	# for (i in 2:dim(cluster.cut)[2]){
		# group_temp<-factor(cluster.cut[,i])
		# group_temp<-as.data.frame(group_temp)
		# adonis_temp<-adonis(dist~group_temp,group_temp,perm=1)
		# Fratio[i]<-adonis_temp$aov.tab$F.Model[1]
	# }
	# return(Fratio)
# }
# data.Fratio<-get.Fratio(data.dist,data.cluster,50)
# plot(data.Fratio,type="h",main="Optimal group number",xlab="group number",ylab="Pseudo F ratio")
validate.cluster<-function(dist,cluster,groupnum,sig){
	distmatrix<-as.matrix(dist)
	cluster.cut<-cutree(cluster,groupnum)
	goodcluster<-T
	for(i in 1:(groupnum-1)){
		indices1<-which(cluster.cut==i)
		for(j in (i+1):groupnum){
			indices2<-which(cluster.cut==j)
			distmatrix_temp<-distmatrix[c(indices1,indices2),c(indices1,indices2)]
			dist_temp<-as.dist(distmatrix_temp)
			group_temp<-cluster.cut[c(indices1,indices2)]
			group_temp<-as.data.frame(group_temp)
			adonis_temp<-adonis(dist_temp~group_temp,group_temp,perm=9999,parallel=clus)
			if(adonis_temp$aov$Pr[1]>sig){
				goodcluster<-F
				break
			}
		}
		if(!goodcluster){
			break
		}
	}
	return(goodcluster)
}
opt.clusternum<-function(dist,cluster,startgroupnum=10,sig=0.01){
	goodcluster<-T
	groupnum<-startgroupnum
	while(goodcluster){
		goodcluster<-validate.cluster(dist,cluster,groupnum,sig)
		groupnum<-groupnum+1
	}
	groupnum<-groupnum-2
	if(groupnum<startgroupnum){
		stop("Too many groups! Please choose a smaller startgroupnum.")
	}
	return(groupnum)
}
data.cluster.cut<-cutree(data.cluster,opt.clusternum(data.dist,data.cluster,10))
n<-attr(data.dist,"Size")
otulabel<-attr(data.dist,"Label")
edge<-matrix(ncol=4,nrow=n*(n-1)/2)
edgecount<-1
for(i in 1:(n-1)){
	for (j in (i+1):n){
		edge[edgecount,1]<-otulabel[i]
		edge[edgecount,2]<-otulabel[j]
		edge[edgecount,3]<-data[i,j]
		edge[edgecount,4]<-pval[i,j]
		edgecount<-edgecount+1
	}
}
stopCluster(clus)
write.table(edge[edge[,4]<0.01,],file="edge.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(data.cluster.cut,file="caggroup.txt",sep="\t",quote=F,col.names=F)