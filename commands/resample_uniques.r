library(doParallel)

n <- 1000
data <- read.table("GD_usearch_pcrstriped_uniques.size")
data<-sort(data[,1],decreasing = T)
cl <- makeCluster(12)
registerDoParallel(cl)
set.seed(1234)
#length(data) is the number of unique sequences, sum(data) is the total sequences
resample <- replicate(n,sample(1:length(data),size = sum(data),replace = T,prob = data))
#transform resample to factor() using table() to count the length(data) unique sequences in sum(data) total sequences
count <- foreach(i=1:n,.combine = 'cbind',.multicombine = T) %dopar% as.data.frame(table(factor(resample[,i],levels=1:length(data))))[,2]
rm(resample)
rmse <- sqrt(rowSums((data-count)^2)/n)
pval <- pt(data/rmse,df=n,lower.tail = F)
write.table(count,"resample.count.txt",quote=F,row.names = F,col.names = F,sep = "\t")
write.table(cbind(data,data/sum(data)*100,data/rmse,pval),file="resample.pval.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(cbind(range(data[pval<0.005]/sum(data)*100),range(data[pval<0.025]/sum(data)*100),range(data[pval<0.05]/sum(data)*100)),file="resample.range.txt",quote=F,sep="\t",row.names=F,col.names=F)
