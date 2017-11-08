#noprecluster resize_moira_filtered.pl mock_moira_pcrstriped.fasta mock_moira.count_table
#noprecluster ~/usearch8 -search_global mock_moira_pcrstriped.resize.fasta -db ~/Illumina_Miseqrun/run04_20141113/run03mock_reference_21.fasta -id 0 -userout mock_moira_pcrstriped.resize.m8 -userfields query+target+id+pairs+gaps+caln+ids+mism+qcov+qrow+trow -strand plus -top_hit_only -maxaccepts 0 -maxrejects 0
#noprecluster sizelengthids.pl mock_moira_pcrstriped.resize.m8
data <- read.table("mock_moira_pcrstriped.resize.m8.sizelengthids",header=F)
1-sum(data[,1]*data[,4])/sum(data[,1]*data[,2])
