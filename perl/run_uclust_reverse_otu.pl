#!/usr/bin/perl
use strict;
use warnings;

my $usage="perl fun_uclust_reverse_otu.pl unique_sorted.fasta\n";
my $unique=shift or die($usage);

mkdir("uclust_reverse_otu");
system("cp $unique uclust_reverse_otu/query.fasta");
chdir("uclust_reverse_otu");

#Make the first (most abundant) seqs in query.fasta become the new target for search to, and add its abundance to size.txt;
open(FILE,"query.fasta");
my $header=<FILE>;
$header=~/size=([0-9]*)/;
open(LOG,">log.txt");
print LOG "SizeOfCentroid\tNumOfOTUs\tNumOfUnmatchedUniques\n";
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
system("parallel_pick_otus_uclust_ref.py -i query.fasta -r target.fasta -o ref/");
system("filter_fasta.py -f query.fasta -o newquery.fasta -s ref/query_failures.txt"); #get unmatched seqs
system("mv newquery.fasta query.fasta"); #use the unmatched seqs to be new query

my $str=`wc -l ref/query_failures.txt`; #count the unmatched seqs
my @numofunmatch=split(/\s/,$str);
print LOG $numofunmatch[0]."\n";

system("rm -r ref/");

while($numofunmatch[0]>0){
	open(FILE,"query.fasta"); #Make the first (most abundant) seqs in query.fasta become the new target for search to, and add its abundance to size.txt;
	my $header=<FILE>;
	$header=~/size=([0-9]*)/;
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
	system("parallel_pick_otus_uclust_ref.py -i query.fasta -r target.fasta -o ref/");
	system("filter_fasta.py -f query.fasta -o newquery.fasta -s ref/query_failures.txt"); #get unmatched seqs
	system("mv newquery.fasta query.fasta"); #use the unmatched seqs to be new query

	my $str=`wc -l ref/query_failures.txt`; #count the unmatched seqs
	@numofunmatch=split(/\s/,$str);
	print LOG $numofunmatch[0]."\n";

	system("rm -r ref/");

}
close(LOG);
close(REP);
