#!/usr/bin/perl
use strict;
use warnings;

my $usage="perl run_threshold.pl summed_threshold.txt\n";
my $file=shift or die($usage);


open(USEARCHOUT,">usearch_threshold_otu_num.txt");
print USEARCHOUT "Uniquemc\tPercent\tNo.ofOTUs\tSeqsremapped\tTp\tTn\tFp\tFn\tMCC\n";
open(QIIMEOUT,">qiime_threshold_otu_num.txt");
print QIIMEOUT "Uniquemc\tPercent\tNo.ofOTUs\tSeqsremapped\tTp\tTn\tFp\tFn\tMCC\n";
open(MOTHUROUT,">mothur_threshold_otu_num.txt");
print MOTHUROUT "Uniquemc\tPercent\tNo.ofOTUs\tSeqsremapped\tTp\tTn\tFp\tFn\tMCC\n";

open(FILE,$file);
while(<FILE>){
	chomp;
	my @line=split(/\t/,$_);
	my $percent=$line[0];
	my $count=$line[1];

	system("mkdir unique_mc$count\_$percent/");
	system("mkdir unique_mc$count\_$percent/usearch");
	system("mkdir unique_mc$count\_$percent/mothur");
    system("cp filtered.resize.fasta filtered.resize.nogap.fasta filtered.resize.count_table unique_mc$count\_$percent");

	chdir("unique_mc$count\_$percent/");

    system("~/usearch8 -sortbysize filtered.resize.nogap.fasta -fastaout filtered.resize.nogap.abundant.fasta -minsize $count");
    system("~/usearch8 -sortbysize filtered.resize.fasta -fastaout filtered.resize.abundant.fasta -minsize $count");
    ###################
    #USEARCH part
    ###################
    system("~/usearch8 -cluster_otus filtered.resize.nogap.abundant.fasta -otus usearch/mock_usearch_pcrstriped_otus.fasta -uparseout usearch/mock_usearch_pcrstriped_otus.log");
    system("fasta_number.py usearch/mock_usearch_pcrstriped_otus.fasta OTU > usearch/mock_usearch_pcrstriped_otus_numbered.fasta");
    system("~/usearch8 -usearch_global filtered.resize.nogap.fasta -db usearch/mock_usearch_pcrstriped_otus_numbered.fasta -strand plus -id 0.97 -uc usearch_readmap.uc");
    ###################
    #QIIME part
    ###################
    system("pick_de_novo_otus.py -i filtered.resize.nogap.abundant.fasta -o qiime -a"); 
    system("~/usearch8 -usearch_global filtered.resize.nogap.fasta -db qiime/rep_set/filtered.resize.nogap.abundant_rep_set.fasta -strand plus -id 0.97 -uc qiime_readmap.uc");
    ##################
    #mothur part
    ##################
    system("mv filtered.resize.abundant.fasta filtered.resize.count_table mothur");
    chdir("mothur");
	system('~/mothur/mothur "#pre.cluster(fasta=filtered.resize.abundant.fasta,count=filtered.resize.count_table,diffs=4,processors=12)"');
	system('~/mothur/mothur "#dist.seqs(fasta=filtered.resize.abundant.precluster.fasta,cutoff=0.2,processors=12)"');
	system('~/mothur/mothur "#cluster(column=filtered.resize.abundant.precluster.dist,count=filtered.resize.abundant.precluster.count_table)"');
	system('~/mothur/mothur "#make.shared(list=filtered.resize.abundant.precluster.an.unique_list.list,count=filtered.resize.abundant.precluster.count_table,label=0.03)"');
	system('~/mothur/mothur "#get.oturep(fasta=filtered.resize.abundant.precluster.fasta,count=filtered.resize.abundant.precluster.count_table,list=filtered.resize.abundant.precluster.an.unique_list.list,column=filtered.resize.abundant.precluster.dist,label=0.03)"');
    system("mv filtered.resize.abundant.precluster.an.unique_list.*.rep.fasta filtered.resize.abundant.precluster.an.unique_list.0.03.rep.fasta");
    system("sed '/^>/!s/-//g' filtered.resize.abundant.precluster.an.unique_list.0.03.rep.fasta >filtered.resize.abundant.precluster.an.unique_list.0.03.rep.sed.fasta");
    system("~/usearch8 -usearch_global ../filtered.resize.nogap.fasta -db filtered.resize.abundant.precluster.an.unique_list.0.03.rep.sed.fasta -strand plus -id 0.97 -uc ../mothur_readmap.uc");
    chdir("../");
    system("rm filtered.resize*");
    chdir("..");
}
close(FILE);

system('ls unique_mc*/*.list |parallel --gnu -j 24 sens_spec.pl {} filtered.resize.dist {.}.sensspec');

open(FILE,$file);
while(<FILE>){
	chomp;
	my @line=split(/\t/,$_);
	my $percent=$line[0];
	my $count=$line[1];
	chdir("unique_mc$count\_$percent/");
    ###################
    #calculate how many OTUs and how many seqs can be remapped to OTUs
    ###################
    my @methodnames=("usearch","qiime","mothur");
    my @remapnum;
    my @otuhashes;
    my @otunumbers;
    foreach my $method(@methodnames){
        my %otuassign;
        open(UC,$method."_readmap.uc");
        my $num=0;
        while(<UC>){
            chomp;
            if(/^H/){
	            my @line=split(/\t/,$_);
	            $line[8]=~ /size=(\d*)/;    #$line[8] is seq name, $line[9] is OTU name
                push(@{$otuassign{$line[9]}},$line[8]);
	            $num+=$1;
            }
        }
        push(@remapnum,$num);
        push(@otuhashes,{%otuassign});
        close(UC);
    ###################
    #make list file
    ###################
        my $listname=$method."_readmap.list";
        open(LIST,">".$listname);
        my $otunum=keys(%otuassign);
        push(@otunumbers,$otunum);
        print LIST "0.03\t$otunum";
        while(my($key,$value)=each(%otuassign)){
            print LIST "\t".join(",",@{$value});
        }
        close(LIST);
    }
    ###################
    #calculate MCC based on the pairwise distance of seqs
    #TP:similarity between seqs>=97% & in the same OTU
    #TN:similarity between seqs<97% & in different OTUs
    #FP:similarity between seqs<97% & in the same OTU
    #FN:similarity between seqs>=97% & in different OTUs
    ###################
    my (@tps,@tns,@fps,@fns,@mccs);
    open(SENSSPEC,"usearch_readmap.sensspec");
    <SENSSPEC>;
    $_=<SENSSPEC>;
    chomp $_;
    @line=split(/\t/,$_);
    push(@tps,$line[0]);
    push(@tns,$line[1]);
    push(@fps,$line[2]);
    push(@fns,$line[3]);
    push(@mccs,$line[4]);
    close(SENSSPEC);
    open(SENSSPEC,"qiime_readmap.sensspec");
    <SENSSPEC>;
    $_=<SENSSPEC>;
    chomp $_;
    @line=split(/\t/,$_);
    push(@tps,$line[0]);
    push(@tns,$line[1]);
    push(@fps,$line[2]);
    push(@fns,$line[3]);
    push(@mccs,$line[4]);
    close(SENSSPEC);
    open(SENSSPEC,"mothur_readmap.sensspec");
    <SENSSPEC>;
    $_=<SENSSPEC>;
    chomp $_;
    @line=split(/\t/,$_);
    push(@tps,$line[0]);
    push(@tns,$line[1]);
    push(@fps,$line[2]);
    push(@fns,$line[3]);
    push(@mccs,$line[4]);
    close(SENSSPEC);
    ###################
    #print the *_threshold_otu_num.txt, including num of OTUs, seqs being remmaped, MCCmock and MCCdist
    ###################
    my $otunum=shift @otunumbers;
    my $remap=shift @remapnum;
    my $tp=shift @tps;
    my $tn=shift @tns;
    my $fp=shift @fps;
    my $fn=shift @fns;
    my $mcc=shift @mccs;
    print USEARCHOUT "$count\t$percent\t$otunum\t$remap\t$tp\t$tn\t$fp\t$fn\t$mcc\n";
    $otunum=shift @otunumbers;
    $remap=shift @remapnum;
    $tp=shift @tps;
    $tn=shift @tns;
    $fp=shift @fps;
    $fn=shift @fns;
    $mcc=shift @mccs;
    print QIIMEOUT "$count\t$percent\t$otunum\t$remap\t$tp\t$tn\t$fp\t$fn\t$mcc\n";
    $otunum=shift @otunumbers;
    $remap=shift @remapnum;
    $tp=shift @tps;
    $tn=shift @tns;
    $fp=shift @fps;
    $fn=shift @fns;
    $mcc=shift @mccs;
    print MOTHUROUT "$count\t$percent\t$otunum\t$remap\t$tp\t$tn\t$fp\t$fn\t$mcc\n";
    chdir("..");
}
