data <- read.table("mock_moira_pcrstriped.resize.m8.sizelengthids")
a <- with(data,sum(V1[V3<80]))
for (i in 81:101){
	a <- c(a,with(data,sum(V1[V3>=(i-1) & V3<i])))
}
a <- rev(a)
write.table(a/sum(a)*100,file="seq_accuracy.txt",row.names=F,col.names=F)
