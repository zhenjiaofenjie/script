#!/usr/bin/perl
use strict;
use warnings;

my $usage="perl convert_mothur_fasta_for_qiime.pl groupsfile fastafile outfile";
my $samplefile=shift;
my $fastafile=shift;
my $outfile=shift;

if(!$samplefile || !$fastafile || !$outfile){
	die($usage);
}

open(SAMPLE,$samplefile);
my %samplehash;
while(<SAMPLE>){
	chomp;
	my @line=split(/\t/);
	$samplehash{$line[0]}=$line[1];
}
close(SAMPLE);

my $qiimenum=10000000;
my $finalseq="";
open(FASTA,$fastafile);
while(<FASTA>){
	chomp;
	my @line=split(/\s+/,$_);
	$line[0]=~ s/^>//;
	if($samplehash{$line[0]}){
		$finalseq.=">".$samplehash{$line[0]}."\_".$qiimenum." ".$line[0]." orig_bc= new_bc= bc_differs=0 quality=\n";
		$_=<FASTA>;
		$finalseq.=$_;
		$qiimenum++;
	}else{
		$_=<FASTA>;
		print $_."\n";
	}
}
close(FASTA);

open(OUT,">".$outfile);
print OUT $finalseq;
close(OUT);
