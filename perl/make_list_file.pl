#!/usr/bin/perl
use warnings;
use strict;

my $outfile=shift;

opendir(DIR,"./");
my %samplefile;
my @dir=readdir(DIR);
foreach my $file(@dir){
	if($file=~ /\.fastq$/){
		my @filename=split(/_/,$file);
		if($file=~ /R1\_001.fastq/){
			${$samplefile{$filename[0]}}[0]=$file;
		}else{
			${$samplefile{$filename[0]}}[1]=$file;
		}
	}
}

open(OUT,">".$outfile);
foreach my $key(keys(%samplefile)){
	print OUT $key."\t".${$samplefile{$key}}[0]."\t".${$samplefile{$key}}[1]."\n";
}
close(OUT);
