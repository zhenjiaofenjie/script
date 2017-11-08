#!/usr/bin/perl
use strict;
use warnings;

my $usage="perl run_AOR_moira.pl threshold\n";
my $count=shift or die($usage);

system("mkdir unique_mc$count/");
system("mkdir unique_mc$count/usearch");
system("mkdir unique_mc$count/mothur");
chdir("unique_mc$count/");
system("~/usearch8 -sortbysize ../filtered.resize.nogap.fasta -fastaout filtered.resize.nogap.abundant.fasta -minsize $count");
system("~/usearch8 -sortbysize ../filtered.resize.fasta -fastaout mothur/filtered.resize.abundant.fasta -minsize $count");
###################
#USEARCH part
###################
system("~/usearch8 -cluster_otus filtered.resize.nogap.abundant.fasta -otus usearch/mock_usearch_pcrstriped_otus.fasta -uparseout usearch/mock_usearch_pcrstriped_otus.log");
system("fasta_number.py usearch/mock_usearch_pcrstriped_otus.fasta OTU > usearch/mock_usearch_pcrstriped_otus_numbered.fasta");
system("~/usearch8 -usearch_global ../filtered.resize.nogap.fasta -db usearch/mock_usearch_pcrstriped_otus_numbered.fasta -strand plus -id 0.97 -uc usearch_readmap.uc");
###################
#QIIME part
###################
system("pick_de_novo_otus.py -i filtered.resize.nogap.abundant.fasta -o qiime -a"); 
system("~/usearch8 -usearch_global ../filtered.resize.nogap.fasta -db qiime/rep_set/filtered.resize.nogap.abundant_rep_set.fasta -strand plus -id 0.97 -uc qiime_readmap.uc");
##################
#mothur part
##################
system("cp ../filtered.resize.count_table mothur");
chdir("mothur");
system('~/mothur/mothur "#pre.cluster(fasta=filtered.resize.abundant.fasta,count=filtered.resize.count_table,diffs=4,processors=12)"');
system('~/mothur/mothur "#dist.seqs(fasta=filtered.resize.abundant.precluster.fasta,cutoff=0.03,processors=12)"');
system('~/mothur/mothur "#cluster(column=filtered.resize.abundant.precluster.dist,count=filtered.resize.abundant.precluster.count_table)"');
system('~/mothur/mothur "#make.shared(list=filtered.resize.abundant.precluster.opti_mcc.list,count=filtered.resize.abundant.precluster.count_table,label=0.03)"');
system('~/mothur/mothur "#get.oturep(fasta=filtered.resize.abundant.precluster.fasta,count=filtered.resize.abundant.precluster.count_table,list=filtered.resize.abundant.precluster.opti_mcc.list,column=filtered.resize.abundant.precluster.dist,label=0.03)"');
system("sed '/^>/!s/-//g' *.rep.fasta >filtered.resize.abundant.precluster.an.unique_list.0.03.rep.sed.fasta");
system("~/usearch8 -usearch_global ../../filtered.resize.nogap.fasta -db filtered.resize.abundant.precluster.an.unique_list.0.03.rep.sed.fasta -strand plus -id 0.97 -uc ../mothur_readmap.uc");
chdir("../");
chdir("..");
