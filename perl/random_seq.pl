#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $USAGE="USAGE: random_seq.pl -fasta/f fastafile -ratio/r ratio_number -out/o number_of_outfile\n";

my ($fasta,$ratio,$outnum,$help)=
   (undef, undef, undef,  "");

&GetOptions('fasta|f=s'	=>\$fasta,
            'ratio|r=f'	=>\$ratio,
            'out|o=i'	=>\$outnum,
            'help|h'	=>\$help,
);

if($help || !$fasta || !$ratio || !$outnum){
	if(!$fasta){
		print "Please verify a fasta file!\n";
	}
	if(!$ratio){
		print "Please verify the ratio number of random picked seqs!\n";
	}
	if(!$outnum){
		print "Please verify how many out files to be generated!\n";
	}
	die($USAGE);
}

my @seqs;
my $seq;
open(FILE,$fasta);
while(<FILE>){
	if(/>/ && $seq){
		push(@seqs,$seq);
		$seq=$_;
	}else{
		$seq=$seq.$_;
	}
}
push(@seqs,$seq);
close(FILE);

$fasta=~ s/\.fasta|\.fa|\.fna//;
my $j=$ratio;
#while($ratio<1){
	my $limit=$ratio*@seqs;
	for(my $i=1;$i<=$outnum;$i++){
		my %hash;
		while((keys %hash)<$limit){
			$hash{int(rand(@seqs))}=1;
		}
		open(OUT,">$fasta\_$j\_$i.fasta");
		foreach(keys %hash){
			print OUT $seqs[$_];
		}
		close(OUT);
	}
	$ratio=$ratio+0.1;
	#$j++;
#}
print "done!\n";
