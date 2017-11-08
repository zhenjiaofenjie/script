n <- 1000
data <- read.table("uniques.size")
data<-sort(data[,1],decreasing = T)
cl <- makeCluster(12)
registerDoParallel(cl)
set.seed(1234)
#length(data) is the number of unique sequences, sum(data) is the total sequences
resample <- replicate(n,sample(1:length(data),size = sum(data),replace = T,prob = data))
#transform resample to factor() using table() to count the length(data) unique sequences in sum(data) total sequences
count <- foreach(i=1:n,.combine = 'cbind',.multicombine = T) %dopar% as.data.frame(table(factor(resample[,i],levels=1:length(data))))[,2]
####################################################
#count <- read.table("resample.count.txt")
###################################################
rm(resample)
r.b <- rowMeans(count)
r.adj <- 2*data-r.b
r.quan <- parApply(cl,count,1,quantile,probs=c(0.005,0.5,0.995),names=F)
r.quan <- t(r.quan)
ci.min <- r.adj-(r.b-r.quan[,1])
ci.max <- r.adj+(r.quan[,3]-r.b)
rmse <- sqrt(rowSums((data-count)^2)/n)
pval <- pt(data/rmse,df=n,lower.tail = F)
########################################################
#out.min <- (data-ci.min)<0 #if real count < ci.min
#out.max <- (data-ci.max)>0 #if real count > ci.max
#out <- out.min+out.max # 1 -> outliers, 0 -> reliable
########################################################
write.table(cbind(data,data/sum(data)*100,r.quan,ci.min,ci.max,pval,data/rmse),"resample.99percentile.txt",quote=F,row.names = F,col.names = c("abund","relative_abund","lower_quan","median","higher_quan","lower_ci","higer_ci","pval","signal/noise"),sep = "\t")
write.table(count,"resample.count.txt",quote=F,row.names = F,col.names = F,sep = "\t")
#write.table(cbind(data,data/sum(data)*100,out),file="resample.out.txt",quote=F,sep="\t",row.names=F,col.names=F)
write.table(cbind(range(data[ci.min<=0]),range(data[ci.min<=0]/sum(data)*100),range(data[ci.min>0]),range(data[ci.min>0]/sum(data)*100)),file="resample.99percentile.range.txt",quote=F,sep="\t",row.names=F,col.names=c("abund_ci<=0","relative_abund_ci<=0","abund_ci>0","relative_abund_ci>0"))
stopCluster(cl)
