#!/usr/bin/perl
use strict;
use warnings;

my @methods=("usearch","qiime","mothur");
foreach my $method(@methods){
##############
#make phylo tree	
###############
	if($method eq "usearch"){
		system("parallel_align_seqs_pynast.py -i $method/mock_usearch_pcrstriped_otus_numbered.fasta -o $method/");
		system("filter_alignment.py -i $method/mock_usearch_pcrstriped_otus_numbered_aligned.fasta -o $method/");
		system("make_phylogeny.py -i $method/mock_usearch_pcrstriped_otus_numbered_aligned_pfiltered.fasta -o $method/rep_set.tre");
	}elsif($method eq "mothur"){
		system("sed -e '/^>/!s/-//g' -e 's/.*\\(Otu\\w*\\).*/>\\1/' $method/filtered.resize.precluster.an.unique_list.0.03.rep.fasta >$method/rep.sed.fasta");
		system("parallel_align_seqs_pynast.py -i $method/rep.sed.fasta -o $method/");
		system("filter_alignment.py -i $method/rep.sed_aligned.fasta -o $method/");
		system("make_phylogeny.py -i $method/rep.sed_aligned_pfiltered.fasta -o $method/rep_set.tre");
		system("cp $method/filtered.resize.precluster.an.unique_list.0.03.biom $method/otu_table.biom");
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
