size<-read.table("GD_usearch_pcrstriped_uniques.size")
xaxis<-sort(unique(size[,1]))/sum(size[,1])*100
summed <- numeric()
j=1
for (i in sort(unique(size[,1]))){
  summed[j]<-sum(size[,1]==i)*i
  j=j+1
}
relative_abund <- summed/sum(size[,1])*100
options(scipen=0)
pdf("summed_size.pdf")
par(mfrow=c(2,1))
plot(xaxis,relative_abund,log="x",xlab="relative abundance (%)",ylab="summed abundance (%)",type="l",xaxt="n")
thresholds<-numeric()
j<-1
for (i in seq(floor(cumsum(relative_abund)[1]),40,1)){
  point<-which(cumsum(relative_abund)>=i & cumsum(relative_abund)<i+1)[1]
  point
  if(!is.na(point)){
#    abline(v=xaxis[point],lty="dashed")
    axis(1,at=xaxis[point],labels=signif(xaxis[point],1),las=2)
    thresholds[j]<-xaxis[point]*sum(size[,1])/100
    names(thresholds)[j]<-paste(i,"%",sep="")
    j<-j+1
  }else{
    next
  }
}
plot(xaxis,cumsum(relative_abund),log="x",xlab="relative abundance (%)",ylab="cumulative abundance (%)",type="l",yaxt="n")
abline(h=c(relative_abund[1],40),lty="dashed")
axis(2,at=c(relative_abund[1],40),labels=c(round(relative_abund[1],1),40))
dev.off()
write.table(thresholds,file="summed_thresholds.txt",col.names = F,quote=F,sep="\t")
