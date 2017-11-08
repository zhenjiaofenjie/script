#!/usr/bin/perl
use strict;
use warnings;

######################
#Change to use Usearch8.0.1623.
#Change maxee from 0.5 to 1.
#-Modified by Jing Wang 20150505.
######################

my $USAGE="perl run_usearch_mergefilter.pl samplefile\n";
my $samplefile=shift or die($USAGE);
my $projectname=$samplefile;
$projectname=~ s/\.txt//;

open(FILE,$samplefile);
while(<FILE>){
	chomp;
	my @line=split(/\t/);
	my $samplename=shift(@line);
	my $fastqr1=shift(@line);
	my $fastqr2=shift(@line);
	system("~/usearch8 -fastq_mergepairs $fastqr1 -reverse $fastqr2 -fastq_minovlen 50 -fastq_truncqual 2 -fastq_minmergelen 400 -fastqout $samplename\_merged.fastq");
	system("~/usearch7.0.1090_i86linux32 -fastq_filter $samplename\_merged.fastq -fastaout $samplename\_filtered.fasta -relabel $samplename\_1000000 -eeout -fastq_maxee 1.0");
 	open(SEDIN,"$samplename\_filtered.fasta");
	open(SEDOUT,">$samplename\_filtered_sed.fasta");
	while(<SEDIN>){
		chomp;
		if(/>/){	
			$_=~ s/$/barcodelabel=$samplename;/;
		}
		print SEDOUT $_."\n";
	}
	close(SEDIN);
	close(SEDOUT);
}
