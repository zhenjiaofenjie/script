#!/usr/bin/perl
use warnings;
use strict;

my $usage="perl sizelengthids.pl uniques_hit.uc
uniqueseqs.uc is the outcome using 
uclust --input uniques.fasta --lib reference.fasta --uc uniques.uc --id 0.00 --libonly
grep \"^H\" uniques.uc >uniques_hit.uc\n";

my $file=shift or die $usage;
my $out=$file.".sizelengthids";
open(FILE,$file);
open(OUT,">$out");
while(<FILE>){
	chomp;
	my @cols=split(/\s/);
	#print $cols[8]."\n";
	$cols[8]=~ /size=([0-9]*)/;
	my $size=$1;
	my $ids=$cols[3];
	my $length=$cols[2];
	print OUT "$size\t$length\t$ids\n";
}
close(FILE);
close(OUT);
