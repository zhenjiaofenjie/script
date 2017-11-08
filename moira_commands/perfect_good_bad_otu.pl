#!/usr/bin/perl
use strict;
use warnings;

open(UNIQUE,"filtered.resize.align.report");
my %uniquehash;
<UNIQUE>;
while(<UNIQUE>){
	chomp;
	my @line=split(/\t/,$_);
	$line[0]=~ /size=(\d*)/;
	$uniquehash{$line[0]}=$1."\t".$line[15];
}
close(UNIQUE);

open(REPORT,"usearch/usearch.align.report");
my %usearchotuhash;
<REPORT>;
while(<REPORT>){
	chomp;
	my @line=split(/\t/,$_);
	$usearchotuhash{$line[0]}=$line[15];
}
close(REPORT);

open(USEARCH,">usearch_otu_unique.txt");
print USEARCH "uniqueid\tsize\tids\totuids\n";
open(READMAP,"usearch/mock_usearch_pcrstriped_readmap.uc");
while(<READMAP>){
	if(/^H/){
		chomp;
		my @line=split(/\t/,$_);
		print USEARCH $line[8]."\t".$uniquehash{$line[8]}."\t".$usearchotuhash{$line[9]}."\n";
	}
}
close(READMAP);
close(USEARCH);

open(REPORT,"qiime/qiime.align.report");
my %qiimeotuhash;
<REPORT>;
while(<REPORT>){
	chomp;
	my @line=split(/\t/,$_);
	$qiimeotuhash{$line[0]}=$line[15];
}
close(REPORT);

open(QIIME,">qiime_otu_unique.txt");
print QIIME "uniqueid\tsize\tids\totuids\n";
open(OTU,"qiime/uclust_picked_otus/filtered.resize.nogap_otus.txt");
while(<OTU>){
	chomp;
	my @line=split(/\t/,$_);
	my $otuname=shift @line;
	while(my $uniquename=shift @line){
		print QIIME $uniquename."\t".$uniquehash{$uniquename}."\t".$qiimeotuhash{$otuname}."\n";
	}
}
close(OTU);
close(QIIME);

open(REPORT,"mothur/mothur.align.report");
my %mothurotuhash;
<REPORT>;
while(<REPORT>){
	chomp;
	my @line=split(/\t/,$_);
	$mothurotuhash{$line[0]}=$line[15];
}
close(REPORT);

open(MOTHUR,">mothur_otu_unique.txt");
print MOTHUR "uniqueid\tsize\tids\totuids\n";
open(LIST,"mothur/filtered.resize.precluster.an.unique_list.list");
while(<LIST>){
	chomp;
	if(/^0.03/){
		my @line=split(/\t/,$_);
		splice(@line,0,2);
		while(my $uniques=shift @line){
			my $otuname;
			my @unique=split(/,/,$uniques);
			foreach my $uniquename(@unique){
				if($mothurotuhash{$uniquename}){
					$otuname=$uniquename;
				}
			}
			foreach my $uniquename(@unique){
				print MOTHUR $uniquename."\t".$uniquehash{$uniquename}."\t".$mothurotuhash{$otuname}."\n";
			}
		}
	}
}
close(LIST);
close(MOTHUR);

