#!/usr/bin/perl
use strict;
use warnings;

my $usage="Usage: get_seqs_in_qiime_otu.pl otuid.txt uclust.txt outseqid\notuid.txt contains otuid in the first column\n";
my $otuidfile=shift or die($usage);
my $uclust=shift or die($usage);
my $output=shift or die($usage);

open(OTU,$otuidfile);
my %otuid;
while(<OTU>){
	chomp;
	my @line=split(/\t/,$_);
	$otuid{$line[0]}=1;
}
close(OTU);

open(OUT,">".$output);
open(UCLUST,$uclust);
while(<UCLUST>){
	chomp;
	my @line=split(/\t/,$_);
	if($otuid{$line[0]}){
		for(my $i=1;$i<@line;$i++){
			print OUT $line[$i]."\n";
		}
	}
}
close(UCLUST);
close(OTU);
