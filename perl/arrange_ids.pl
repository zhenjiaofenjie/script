#!/usr/bin/perl -w

my $sample = shift or die("USAGE: arrange_dis.pl sample.txt rawid.txt\n");

my $rawid = shift or die("USAGE: arrange_dis.pl sample.txt rawid.txt\n");

open SAM,$sample or die $0;

open RAW,$rawid or die $0;

my %sample_idx;

foreach (<SAM>){
	chomp;
	my $line = $_;
	my @arr = split("\t",$line);
	$sample_idx{$arr[0]} = $arr[1];
	$sample_idx{$arr[0]} -> {TAG} = $arr[2];
}

close SAM;

my %raw_ids;

foreach (<RAW>){
	chomp;
	my $line = $_;
	my @arr = split("\t",$line);
	print $sample_idx{$arr[0]}."\t".$arr[1]."\t".$sample_idx{$arr[0]}->{TAG}."\n" if $sample_idx{$arr[0]};
}

close RAW;