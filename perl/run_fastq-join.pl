#!/usr/bin/perl
use strict;
use warnings;

my $USAGE="perl run_fastq-join.pl samplefile\n";
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
	#system("fastq-join $fastqr1 $fastqr2 -o $samplename");
	system("fastq-join -m 50 $fastqr1 $fastqr2 -o $samplename");
}
