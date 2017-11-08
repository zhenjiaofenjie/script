#!/usr/bin/perl
use warnings;
use strict;

if(!@ARGV){
	die("USAGE: FakeFastq.pl seqfile.fasta\n");
}
my $fasta = $ARGV[0];
my ($fast_head,$fastq,$bases);
if($fasta =~ /(\S+)\.[fasta fna fa]/){
	$fast_head = $1;
	$fastq = $fast_head.".fastq";
}
else {print "Only files with .fasta can be processed.\n";}

open FASTA,$fasta;
open FASTQ,">$fastq";
select FASTQ;
foreach (<FASTA>){
	chomp;
	if(/>/){
		s/>/\@/;
		print $_."\n";
	}
	else {
		print $_."\n";
		$bases = length($_);
		print "+\n";
		for (1..$bases){
			print "I";
		}
		print "\n";
	}
}
close FASTA;
close FASTQ;
