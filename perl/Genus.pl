#!/usr/bin/perl
use strict;
use 5.010;

my $hierarchy = $ARGV[0] or die("USAGE: Genus.pl rdpfile.txt\n");

my $hieout = $hierarchy;
$hieout =~ s/\.txt/_genus\.txt/;

open HIE,$hierarchy or die("cann't open hierarchy file\n");
open OUT,">$hieout";
select OUT;

foreach (<HIE>){
	chomp;
	if(/genus \"?(?<Genus>[\w\/]+)\"? \((?<Num>\d+)\)/){
		print $+{Genus}."	".$+{Num}."\n";
	}
}
close HIE;
close OUT;