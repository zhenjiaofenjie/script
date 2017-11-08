#totalcount <- 254086
#totalcount <- 130940
#totalcount <- 190818
usearch <- read.table("usearch_otu_summary.txt",header=T)
denovo <- read.table("denovo_otu_summary.txt",header=T)
chimeraslayer <- read.table("denovo_nochimeric_otu_summary.txt",header=T)
uchime <- read.table("denovo_uchime_nochimeric_otu_summary.txt",header=T)
ii <- order(usearch$X.discard)
usearch <- usearch[ii,]
denovo <- denovo[ii,]
chimeraslayer <- chimeraslayer[ii,]
uchime <- uchime[ii,]
#######################
totalcount <- denovo$No.ofseqs[1]
pdf("otu_threshold.pdf")
par(mfcol=c(2,1))
######################
plot(denovo$MinCount/totalcount*100,denovo$No.ofOTU,type="l",col="red",lwd=2,xlab="relative abundance (%)",ylab = "No. of OTUs",log="x")
#######################
#title(main="PWS")
#######################
lines(uchime$MinCount/totalcount*100,uchime$No.ofOTU,type="l",col="purple",lwd=2)
lines(usearch$MinCount/totalcount*100,usearch$No.ofOTU,type="l",col="blue",lwd=2)
lines(chimeraslayer$MinCount/totalcount*100,chimeraslayer$No.ofOTU,type="l",col="green",lwd=2)
legend("topright",c("Usearch","Qiime","Qiime+ChimeraSlayer","Qiime+Uchime"),col=c("blue","red","green","purple"),lwd=2,cex=0.8)
plot(denovo$MinCount/totalcount*100,denovo$No.ofseqs/totalcount*100,type="l",col="red",lwd=2,xlab="relative abundance (%)",ylab = "remapped sequences (%)",ylim=c(50,100),log="x")
title(main="% of sequences set aside",font.main=1,cex.main=1,line=2.5)
axis(3,at=denovo$MinCount/totalcount*100,labels=F)
axis(3,at=denovo$MinCount[c(1,which(denovo$X.discard %in% c(20,25,30,35)),length(denovo$MinCount))]/totalcount*100,labels=denovo$X.discard[c(1,which(denovo$X.discard %in% c(20,25,30,35)),length(denovo$X.discard))])
lines(uchime$MinCount/totalcount*100,uchime$No.ofseqs/totalcount*100,type="l",col="purple",lwd=2)
lines(usearch$MinCount/totalcount*100,usearch$No.ofseqs/totalcount*100,type="l",col="blue",lwd=2)
lines(chimeraslayer$MinCount/totalcount*100,chimeraslayer$No.ofseqs/totalcount*100,type="l",col="green",lwd=2)
#legend("bottomright",c("Usearch","Qiime","Qiime+ChimeraSlayer","Qiime+Uchime"),col=c("blue","red","green","purple"),lwd=2,cex=0.8)
dev.off()
