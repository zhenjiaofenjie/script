#!/usr/bin/perl
use strict;
use warnings;

my $usage="perl fun_usearch_reverse_otu.pl unique_sorted.fasta\n";
my $unique=shift or die($usage);

mkdir("usearch_reverse_otu");
system("cp $unique usearch_reverse_otu/query.fasta");
chdir("usearch_reverse_otu");

#Make the first (most abundant) seqs in query.fasta become the new target for search to, and add its abundance to size.txt;
open(FILE,"query.fasta");
my $header=<FILE>;
$header=~/size=([0-9]*)/;
open(LOG,">log.txt");
print LOG "SizeOfCentroid\tNumOfOTUs\tNumOfUnmatchedUniques\tNumOfMatchedUniques\n";
print LOG $1."\t";
my $seq=<FILE>;
while(<FILE>){
	if($_=~/^>/){
		last;
	}
	$seq.=$_;
}
my $otucount=1;
print LOG $otucount."\t";
open(TAR,">target.fasta");
open(REP,">rep_set.fasta");
print TAR ">OTU$otucount\n".$seq;
print REP ">OTU$otucount\n".$seq;
close(TAR);
close(FILE);

#Search query.fasta on the new target.fasta
system("~/usearch8 -usearch_global query.fasta -db target.fasta -strand plus -id 0.97 -uc readmap.uc");
system(`awk -F "\t" '(/^N/){print \$9}' readmap.uc >newquery.txt`); #get the unmatched seqs
system(`awk -F "\t" '(/^H/){print}' readmap.uc >final_readmap.txt`); #get the matched seqs
system("filter_fasta.py -f query.fasta -o newquery.fasta -s newquery.txt");
system("mv newquery.fasta query.fasta"); #use the unmatched seqs to be new query

my $str=`wc -l newquery.txt`; #count the unmatched seqs
my @numofunmatch=split(/\s/,$str);
print LOG $numofunmatch[0]."\t";

my $numofmatch=`awk -F "\t" 'BEGIN{count=0}(/^H/){count=count+1}END{print count}' readmap.uc`;
print LOG $numofmatch."\n";

while($numofunmatch[0]>0){
	open(FILE,"query.fasta"); #Make the first (most abundant) seqs in query.fasta become the new target for search to, and add its abundance to size.txt;
	my $header=<FILE>;
	$header=~/size=([0-9]*)/;
#	if($1==1){
#		last;  #Avoid singleton to be centroids, due to the time cost.
#	}
	print LOG $1."\t";
	my $seq=<FILE>;
	while(<FILE>){
		if($_=~/^>/){
			last;
		}
		$seq.=$_;
	}
	$otucount++;
	print LOG $otucount."\t";
	open(TAR,">target.fasta"); #target.fasta only contains one sequence
	print TAR ">OTU$otucount\n".$seq;
	print REP ">OTU$otucount\n".$seq; #rep_set.fasta contains all the centroids so far
	close(TAR);
	close(FILE);

	system("~/usearch8 -usearch_global query.fasta -db target.fasta -strand plus -id 0.97 -uc readmap.uc"); #Search query.fasta on the new target.fasta
	system(`awk -F "\t" '(/^N/){print \$9}' readmap.uc >newquery.txt`); #get the unmatched seqs
	system(`awk -F "\t" '(/^H/){print}' readmap.uc >>final_readmap.txt`); #get the matched seqs
	system("filter_fasta.py -f query.fasta -o newquery.fasta -s newquery.txt");
	system("mv newquery.fasta query.fasta"); #use the unmatched seqs to be new query

	my $str=`wc -l newquery.txt`; #count the unmatched seqs
	@numofunmatch=split(/\s/,$str);
	print LOG $numofunmatch[0]."\t";
	
	my $numofmatch=`awk -F "\t" 'BEGIN{count=0}(/^H/){count=count+1}END{print count}' readmap.uc`;
	print LOG $numofmatch."\n";
}
#print LOG "1\t".$otucount+$numofunmatch[0]."\t0";
close(LOG);
close(REP);
