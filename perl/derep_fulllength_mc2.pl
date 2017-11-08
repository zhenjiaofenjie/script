#!/usr/bin/perl
use strict;
use warnings;

my $usage="perl derep_fulllength.pl fastain\n";
my $fastain=shift or die($usage);
my $fastaout=$fastain;
$fastaout=~ s/(\.fasta|\.fa|\.fna)/_uniques_mc2.fasta/;
my $uniquemap=$fastain;
$uniquemap=~ s/(\.fasta|\.fa|\.fna)/_unimap_mc2.txt/;

open(FASTA,$fastain);
my $seqname=<FASTA>;
chomp $seqname;
$seqname=substr($seqname,1);
my $seq="";
my %uniquehash;   #@{$uniquehash{$seq}}=(uniquename,size,seqname in this uniqseq);
while(<FASTA>){
	chomp;
	if(/^>/){
		if(!$uniquehash{$seq}){
			${$uniquehash{$seq}}[0]=$seqname;
			${$uniquehash{$seq}}[2]=$seqname;
		}else{
			${$uniquehash{$seq}}[2].="\t".$seqname;
		}
		${$uniquehash{$seq}}[1]++;
		$seqname=$_;
		$seqname=substr($seqname,1);
		$seq="";
	}else{
		$seq.=$_;
	}
}
if(!$uniquehash{$seq}){
	${$uniquehash{$seq}}[0]=$seqname;
	${$uniquehash{$seq}}[2]=$seqname;
}else{
	${$uniquehash{$seq}}[2].="\t".$seqname;
}
${$uniquehash{$seq}}[1]++;
close(FASTA);

open(OUT,">".$fastaout);
open(MAP,">".$uniquemap);
foreach my $seq(keys %uniquehash){
	if(${$uniquehash{$seq}}[1]>1){
		print OUT ">".${$uniquehash{$seq}}[0]."size=".${$uniquehash{$seq}}[1].";\n".$seq."\n";
		print MAP ${$uniquehash{$seq}}[0]."size=".${$uniquehash{$seq}}[1].";\t".${$uniquehash{$seq}}[2]."\n";
	}
}
close(OUT);
close(MAP);
