#!/usr/bin/perl
use strict;
use warnings;

###################
#read count table of unique seqs
###################
my %seqcount;
open(COUNT,"../filtered.resize.count_table");
while(<COUNT>){
	chomp;
	my @line=split(/\t/,$_);
	my $seqname=shift @line;
	shift @line;
	@{$seqcount{$seqname}}=@line;
}
close(COUNT);

my @methods=("usearch","qiime","mothur");
foreach my $method(@methods){
####################
#reading * readmap.uc
####################
	open(UC,$method."_readmap.uc");
	my %otuassign;
	while(<UC>){
		chomp;
		if(/^H/){
			my @line=split(/\t/,$_);
			push(@{$otuassign{$line[9]}},$line[8]);
		}
	}
	close(UC);

####################
#make OTU table
####################
	open(OTU,">".$method."/otu_table.txt");
	print OTU "Representative_Sequence\t".join("\t",@{$seqcount{'Representative_Sequence'}})."\n";
	while((my $key,my $value)=each(%otuassign)){
		my @abund;
		foreach my $seqname(@{$value}){
			for(my $i=0;$i<@{$seqcount{$seqname}};$i++){
				$abund[$i]+=${$seqcount{$seqname}}[$i];
			}
		}
		print OTU $key."\t".join("\t",@abund)."\n";
	}
	close(OTU);
	system("biom convert -i $method/otu_table.txt -o $method/otu_table.biom --to-json");
###############
#make phylo tree	
###############
	if($method eq "usearch"){
		system("parallel_align_seqs_pynast.py -i $method/mock_usearch_pcrstriped_otus_numbered.fasta -o $method/");
		system("filter_alignment.py -i $method/mock_usearch_pcrstriped_otus_numbered_aligned.fasta -o $method/");
		system("make_phylogeny.py -i $method/mock_usearch_pcrstriped_otus_numbered_aligned_pfiltered.fasta -o $method/rep_set.tre");
	}elsif($method eq "mothur"){
		system("parallel_align_seqs_pynast.py -i mothur/filtered.resize.abundant.precluster.an.unique_list.0.03.rep.sed.fasta -o $method/");
		system("filter_alignment.py -i $method/filtered.resize.abundant.precluster.an.unique_list.0.03.rep.sed_aligned.fasta -o $method/");
		system("make_phylogeny.py -i $method/filtered.resize.abundant.precluster.an.unique_list.0.03.rep.sed_aligned_pfiltered.fasta -o $method/rep_set.tre");
	}
###############
#calc alpha diversity	
###############
	system("alpha_diversity.py -i $method/otu_table.biom -o $method/alpha.txt -m observed_otus,chao1,shannon,simpson");	
###############
#calc beta diversity	
###############
	system("beta_diversity.py -i $method/otu_table.biom -o $method/beta -m euclidean,bray_curtis,unweighted_unifrac,weighted_unifrac -t $method/rep_set.tre");
}
