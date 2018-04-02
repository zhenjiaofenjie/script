library(doParallel)
#Calculate whether the 99% confidential interval of unique sequence contains 0. If so, the unique sequence is not statistically reliable.
n <- 1000 #Set the No. of replicates in bootstrap process.
data <- read.table("uniques.size") #uniques.size is a file contains only one column, in which are the absolute abundances of each unique sequences. 
data<-sort(data[,1],decreasing = T)
cl <- makeCluster(12) #Set the No. of proccessors used during parallel calculation.
registerDoParallel(cl)
set.seed(1234) #The seed of randomized resampling.
#length(data) is the number of unique sequences, sum(data) is the total sequences
resample <- replicate(n,sample(1:length(data),size = sum(data),replace = T,prob = data))
#transform resample to factor() using table() to count the length(data) unique sequences in sum(data) total sequences
count <- foreach(i=1:n,.combine = 'cbind',.multicombine = T) %dopar% as.data.frame(table(factor(resample[,i],levels=1:length(data))))[,2]
rm(resample)
r.b <- rowMeans(count) #Mean of bootstrapped abundance
r.adj <- 2*data-r.b #Adjusted mean of bootstrapped abundance
r.quan <- parApply(cl,count,1,quantile,probs=c(0.005,0.5,0.995),names=F) #99% percentile of bootstrapped abundance
r.quan <- t(r.quan)
ci.min <- r.adj-(r.b-r.quan[,1]) #99% confidantial interval of abundance
ci.max <- r.adj+(r.quan[,3]-r.b)
rmse <- sqrt(rowSums((data-count)^2)/n)
pval <- pt(data/rmse,df=n,lower.tail = F)
#The result of each unique sequence
write.table(cbind(data,data/sum(data)*100,r.quan,ci.min,ci.max,pval,data/rmse),"resample.99percentile.txt",quote=F,row.names = F,col.names = c("abund","relative_abund","lower_quan","median","higher_quan","lower_ci","higer_ci","pval","signal/noise"),sep = "\t")
#The bootstrapped abundance
write.table(count,"resample.count.txt",quote=F,row.names = F,col.names = F,sep = "\t")
#Determine the abundance range of reliable and unreliable sequences.
write.table(cbind(range(data[ci.min<=0]),range(data[ci.min<=0]/sum(data)*100),range(data[ci.min>0]),range(data[ci.min>0]/sum(data)*100)),file="resample.99percentile.range.txt",quote=F,sep="\t",row.names=F,col.names=c("abund_ci<=0","relative_abund_ci<=0","abund_ci>0","relative_abund_ci>0"))
stopCluster(cl)
