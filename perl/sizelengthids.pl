#!/usr/bin/perl
use warnings;
use strict;

my $usage="perl sizelengthids.pl uniqueseqs.m8
uniqueseqs.m8 is the outcome using 
~/usearch8 -search_global inputfile -db reference -id 0 -userout usearchout.m8 -userfields query+target+id+pairs+gaps+caln+ids+mism+qcov+qrow+trow -strand plus -top_hit_only -maxaccepts 0 -maxrejects 0\n";

my $file=shift or die $usage;
my $out=$file.".sizelengthids";
open(FILE,$file);
open(OUT,">$out");
while(<FILE>){
	chomp;
	my @cols=split(/\t/);
	$cols[0]=~ /size=([0-9]*)/;
	my $size=$1;
	my $id=$cols[2];
	my $length=$cols[3];
    my $ids=$cols[6];
	print OUT "$size\t$length\t$id\t$ids\n";
}
close(FILE);
close(OUT);
