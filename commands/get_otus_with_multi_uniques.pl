#!/usr/bin/perl
use strict;
use warnings;

my $usage="get_otus_with_multi_uniques.pl readmap.uc\n";
my $readmap=shift or die($usage);
my $out="otuid_with_multi_uniques.txt";

my %otuid;
open(READMAP,$readmap);
while(<READMAP>){
	chomp;
	my @col=split(/\t/,$_);
	if($otuid{$col[9]}){
		$otuid{$col[9]}++;
	}else{
		$otuid{$col[9]}=1;
	}
}
close(READMAP);

open(OUT,">".$out);
while((my $key, my $value)=each(%otuid)){
	if($value>1){
		print OUT $key."\t".$value."\n";
	}
}
close(OUT)
