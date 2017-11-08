#!/usr/bin/perl
use strict;
use warnings;

######################
#
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
	my $bam=shift(@line);
 	system("samtools bam2fq $bam >$samplename\.fastq");
	system("~/usearch7.0.1090_i86linux32 -fastq_filter $samplename\.fastq -relabel $samplename\_ -eeout -fastaout $samplename\_filtered.fasta -fastq_maxee 4 -fastq_trunclen 150");
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
