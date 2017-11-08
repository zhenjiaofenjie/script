#!/usr/bin/perl
use warnings;
use strict;

my $usage="get_uniques_sizeids_mothur.pl unique.m8 unique.count_table unique2otu.txt rep_set.m8\n";
my $idfile=shift or die($usage);
my $countfile=shift or die($usage);
my $unique2otufile=shift or die($usage);
my $repsetidfile=shift or die($usage);

my %id;
open(ID,$idfile);
while(<ID>){
	chomp;
	my @line=split(/\t/,$_);
	$id{$line[0]}=$line[2];
}
close(ID);

open(COUNT,$countfile);
while(<COUNT>){
	chomp;
	my @line=split(/\t/,$_);
	if($id{$line[0]}){
		$id{$line[0]}.="\t".$line[1];
	}
}
close(COUNT);

open(UNI2OTU,$unique2otufile);
my %unique2otu;
while(<UNI2OTU>){
	chomp;
	my @line=split(/\t/,$_);
	$unique2otu{$line[0]}=$line[1];
}
close(UNI2OTU);

open(REPID,$repsetidfile);
my %repid;
while(<REPID>){
	chomp;
	my @line=split(/\t/,$_);
	$repid{$line[0]}=$line[2];
}
close(REPID);

open(PERFECT,">perfect.sizeids");
open(GOOD,">good.sizeids");
open(BAD,">bad.sizeids");
while(my($key,$value)=each(%id)){
	my $otuid=$repid{$unique2otu{$key}};
	if($otuid==100){
		print PERFECT $value."\n";
	}elsif($otuid>=97){
		print GOOD $value."\n";
	}else{
		print BAD $value."\n";
	}
}
close(PERFECT);
close(GOOD);
close(BAD);
