#!/usr/bin/perl
use strict;
use warnings;

###################
#Added support for moira.py.
#-20160216
##################

my $usage="perl derep_fulllength.pl fastain\n";
my $fastain=shift or die($usage);
my $fastaout=$fastain;
$fastaout=~ s/(\.fasta|\.fa|\.fna)/_uniques.fasta/;
my $uniquemap=$fastain;
$uniquemap=~ s/(\.fasta|\.fa|\.fna)/_unimap.txt/;

open(FASTA,$fastain);
my $seqname=<FASTA>;
chomp $seqname;
$seqname=substr($seqname,1);
my $seqsize="";
if($seqname=~ /size=(\d*)/){   # In the case of moira.py, the seqs have been dereplicated, thus have size annotation.
    $seqsize=$1;
}else{
    $seqsize="";
}
my $seq="";
my %uniquehash;   #@{$uniquehash{$seq}}=(uniquename,size,seqname in this uniqseq);
while(<FASTA>){
	chomp;
	if(/^>/){
		if(!$uniquehash{$seq}){
			${$uniquehash{$seq}}[2]=$seqname;
		    $seqname=~ s/size=\d*;//;
			${$uniquehash{$seq}}[0]=$seqname;
		}else{
			${$uniquehash{$seq}}[2].="\t".$seqname;
		}
        if($seqsize){
            ${$uniquehash{$seq}}[1]+=$seqsize;
        }else{
		    ${$uniquehash{$seq}}[1]++;
        }
		$seqname=$_;
		$seqname=substr($seqname,1);
		$seq="";
		if($seqname=~ /size=(\d*)/){
		    $seqsize=$1;
		}else{
            $seqsize="";
        }
	}else{
		$seq.=$_;
	}
}
if(!$uniquehash{$seq}){
	${$uniquehash{$seq}}[2]=$seqname;
    $seqname=~ s/size=\d*;//;
	${$uniquehash{$seq}}[0]=$seqname;
}else{
	${$uniquehash{$seq}}[2].="\t".$seqname;
}
if($seqsize){
    ${$uniquehash{$seq}}[1]+=$seqsize;
}else{
    ${$uniquehash{$seq}}[1]++;
}

close(FASTA);

open(OUT,">".$fastaout);
open(MAP,">".$uniquemap);
foreach my $seq(keys %uniquehash){
	print OUT ">".${$uniquehash{$seq}}[0]."size=".${$uniquehash{$seq}}[1].";\n".$seq."\n";
	print MAP ${$uniquehash{$seq}}[0]."size=".${$uniquehash{$seq}}[1].";\t".${$uniquehash{$seq}}[2]."\n";
}
close(OUT);
close(MAP);
