#!/usr/bin/perl
use warnings;
use strict;

my $usage="filter_contaminant.pl contaiminant.txt file out\n";
my $contaminantfile=shift or die($usage);
my $file=shift or die($usage);
my $out=shift or die($usage);

open(CON,$contaminantfile);
my %conhash;
while(<CON>){
	chomp;
	my @line=split(/\t/,$_);
	$conhash{$line[0]}=1;
}
close(CON);

open(FILE,$file);
open(OUT,">$out");
while(<FILE>){
	my @line=split(/\t/,$_);
	if($conhash{$line[0]}){
		next;
	}else{
		print OUT;
	}
}
close(FILE);
close(OUT);
