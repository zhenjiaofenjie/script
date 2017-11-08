#!/usr/bin/perl
use strict;
use warnings;

my $usage="USAGE: perl SeqErrorStat.pl usearchout.m8
usearchout.m8 is the outcome using\n~/usearch8 -search_global inputfile -db reference -id 0 -userout usearchout.m8 -userfields query+target+id+pairs+gaps+caln+ids+mism+qcov+qrow+trow -strand plus -top_hit_only -maxaccepts 0 -maxrejects 0\n";
my $usearchout=shift or die $usage;
my $seqout=$usearchout.".seq.stat";
my $positionout=$usearchout.".position.stat";
my $seqnlineout=$usearchout.".seqnline.stat";
my $summary=$usearchout.".summary.stat";
my $gccontent=$usearchout.".gc.stat";

open(FILE,$usearchout);
open(SEQOUT,">".$seqout);
open(SEQNLINEOUT,">".$seqnlineout);
print SEQOUT "queryid\ttargetid\tA2T\tA2G\tA2C\tA2gap\tT2A\tT2G\tT2C\tT2gap\tC2A\tC2T\tC2G\tC2gap\tG2A\tG2T\tG2C\tG2gap\tgap2A\tgap2T\tgap2G\tgap2C\n";
print SEQNLINEOUT "queryid\ttargetid\tstatus\tposition\n";
my @char=("A","T","G","C","-");
my %position;
my %summary;
my %gc_orig;
my %gc_seq;
while(<FILE>){
	chomp;
	my @line=split(/\t/,$_);
	my %hash;
	print SEQOUT $line[0]."\t".$line[1];
	for(my $i=0;$i<length($line[9]);$i++){
		my $querychar=substr($line[9],$i,1);
		my $targetchar=substr($line[10],$i,1);
		${$gc_orig{$targetchar}}[$i]++;
		${$gc_seq{$querychar}}[$i]++;
		if($querychar ne $targetchar){
			${$hash{$targetchar}}{$querychar}+=1;
			${${$position{$targetchar}}{$querychar}}[$i]+=1;
			print SEQNLINEOUT $line[0]."\t".$line[1]."\t$targetchar->$querychar\t$i\n";
		}
	}
	foreach my $from(@char){
		foreach my $to(@char){
			if ($from ne $to){
				if(${$hash{$from}}{$to}){
					print SEQOUT "\t".${$hash{$from}}{$to};
					${${$summary{$line[1]}}{$from}}{$to}+=${$hash{$from}}{$to};
					${${$summary{"all"}}{$from}}{$to}+=${$hash{$from}}{$to};
					
				}else{
					print SEQOUT "\t0";
				}
			}
		}
	}
	print SEQOUT "\n";
}
close(FILE);
close(SEQOUT);
close(SEQNLINEOUT);

open(SUMMARY,">".$summary);
print SUMMARY "targetid\tA2T\tA2G\tA2C\tA2gap\tT2A\tT2G\tT2C\tT2gap\tC2A\tC2T\tC2G\tC2gap\tG2A\tG2T\tG2C\tG2gap\tgap2A\tgap2T\tgap2G\tgap2C\n";
foreach my $key(sort keys %summary){
	print SUMMARY $key;
	foreach my $from(@char){
		foreach my $to(@char){
			if ($from ne $to){
				if (${${$summary{$key}}{$from}}{$to}){
					print SUMMARY "\t".${${$summary{$key}}{$from}}{$to};
				}else{
					print SUMMARY "\t0";
				}
			}
		}
	}
	print SUMMARY "\n";
}
close(SUMMARY);

open(POSOUT,">".$positionout);
foreach my $from(@char){
	foreach my $to(@char){
		if ($from ne $to){
			foreach my $pos(@{${$position{$from}}{$to}}){if(!$pos){$pos=0;}}
			print POSOUT "$from->$to\t@{${$position{$from}}{$to}}\n";
		}
	}
}
close(POSOUT);

open(GC,">".$gccontent);
print GC "Original\n";
foreach my $key(sort keys(%gc_orig)){
	if($key ne "-"){
		foreach my $pos(@{$gc_orig{$key}}){
			if(!$pos){
				$pos=0;
			}
		}
		print GC "$key\t@{$gc_orig{$key}}\n";
	}
}
print GC "Sequence\n";
foreach my $key(sort keys(%gc_seq)){
        if($key ne "-"){
		foreach my $pos(@{$gc_seq{$key}}){
			if(!$pos){
				$pos=0;
			}
		}
		print GC "$key\t@{$gc_seq{$key}}\n";
	}
}
close(GC);
