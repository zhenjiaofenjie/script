#!/usr/bin/perl
use strict;

my $file = $ARGV[0] or die("USAGE: ClearFa.pl seqfile.fasta\n");

my $outfile = $file;
$outfile =~ s/fasta/fa/;

open File,$file;
open OUT,">$outfile";
select OUT;

foreach(<File>){
	unless(/^\>/){
	s/\.//g;
	s/\-//g;
	print;
}
else {
	print;
}
}
close File;
close OUT;