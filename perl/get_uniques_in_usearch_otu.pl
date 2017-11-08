#!/usr/bin/perl
use strict;
use warnings;

my $usage="Usage: get_uniques_in_usearch_otu.pl otuid.txt readmap.uc unique.m8 outunique.m8\notuid.txt contains otuid in the first column\n";
my $otuidfile=shift or die($usage);
my $readmapfile=shift or die($usage);
my $uniquem8=shift or die($usage);
my $output=shift or die($usage);

open(OTU,$otuidfile);
my %otuid;
while(<OTU>){
	chomp;
	my @line=split(/\t/,$_);
	$otuid{$line[0]}=1;
}
close(OTU);

open(READMAP,$readmapfile);
my %uniqueid;
while(<READMAP>){
	chomp;
	my @line=split(/\t/,$_);
	if($line[0] eq "H" && $otuid{$line[9]}){
		$uniqueid{$line[8]}=$line[9];
	}
}
close(READMAP);

open(OUT,">".$output);
open(UNIQUEM8,$uniquem8);
while(<UNIQUEM8>){
	chomp;
	my @line=split(/\t/,$_);
	if($uniqueid{$line[0]}){
		print OUT $_."\t".$uniqueid{$line[0]}."\n";
	}
}
close(UNIQUEM8);
close(OUT);
