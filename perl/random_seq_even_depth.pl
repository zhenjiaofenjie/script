#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $USAGE="USAGE: random_seq.pl -fasta/f fastafile -even/e number_of_seqs_per_sample -out/o number_of_outfile\n";

my ($fasta,$even,$outnum,$help)=
   (undef, undef, undef,  "");

&GetOptions('fasta|f=s'	=>\$fasta,
            'even|e=i'	=>\$even,
            'out|o=i'	=>\$outnum,
            'help|h'	=>\$help,
);

if($help || !$fasta || !$even || !$outnum){
	if(!$fasta){
		print "Please verify a fasta file!\n";
	}
	if(!$even){
		print "Please verify the even number of random picked seqs!\n";
	}
	if(!$outnum){
		print "Please verify how many out files to be generated!\n";
	}
	die($USAGE);
}

my %seqs;
my $seq;
my $samplename;
open(FILE,$fasta);
while(<FILE>){
	if(/>/){
		if($seq){
			push(@{$seqs{$samplename}},$seq);
		}
		$seq=$_;
		$_=~ /^>(.*)_/;
		$samplename=$1;
	}else{
		$seq=$seq.$_;
	}
}
push(@{$seqs{$samplename}},$seq);
close(FILE);

$fasta=~ s/\.fasta|\.fa|\.fna//;
for(my $i=1;$i<=$outnum;$i++){
	open(OUT,">$fasta\_even$even\_$i.fasta");
	foreach my $samplename(keys(%seqs)){
		my %hash;
		my $seqnum=@{$seqs{$samplename}};
		if($seqnum<$even){
			next;
		}elsif($seqnum==$even){
			print OUT join("",@{$seqs{$samplename}});
			next;
		}
		while((keys %hash)<$even){
			$hash{int(rand($seqnum))}=1;
		}
		foreach(keys %hash){
			print OUT ${$seqs{$samplename}}[$_];
		}
	}
	close(OUT);
}
print "done!\n";
