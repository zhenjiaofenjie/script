data <- read.table("mothur.nocontaminant.align.report",header=T)
a <- hist(data$SimBtwnQuery.Template,plot=F,breaks=c(0,90:101),right=F)
write.table(a$counts,file="otu_best_hit.txt",col.names=F,row.names=F)
