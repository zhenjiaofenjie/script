#!/usr/bin/perl
use strict;
use warnings;

my $usage="Usage: perl alpha_split.pl alpha_paste.txt\n";
my $alpha=shift or die($usage);

open(ALPHA,$alpha);
open(OS,">observed_species.txt");
open(CHAO1,">chao1.txt");
open(SHANNON,">shannon.txt");
open(SIMPSON,">simpson.txt");

while(<ALPHA>){
	chomp;
	my @line=split(/\t/,$_);
	for (my $i=0;$i<9;$i++){
		print OS $line[5*$i+1]."\t";
		print CHAO1 $line[5*$i+2]."\t";
		print SHANNON $line[5*$i+3]."\t";
		print SIMPSON $line[5*$i+4]."\t";
	}
	print OS $line[46]."\n";
        print CHAO1 $line[47]."\n";
        print SHANNON $line[48]."\n";
        print SIMPSON $line[49]."\n";
}
close(ALPHA);
close(OS);
close(CHAO1);
close(SHANNON);
close(SIMPSON);
