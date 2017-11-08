#!/usr/bin/perl
use strict;
use warnings;

my $usage="sens_spec.pl listfile distfile output\n";
my $list=shift or die($usage);
my $dist=shift or die($usage);
my $out=shift or die($usage);

my($Tp,$Tn,$Fp,$Fn,$MCC)=
(0, 0,   0,  0,  0);

print "Reading list file...\n";
open(LIST,$list);
my %otuassign;
while(<LIST>){ #first suppose all pairs of seqs has similarity<97%
	if(/^0.03/){
		chomp;
		my @seqnames=split(/\t|,/,$_); #count seqs
		my $seqcount=@seqnames-2; #acrual number of seqs
		@seqnames=();
		my @line=split(/\t/,$_);
		splice(@line,0,2);
		my $otuid=1;
		foreach my $otu(@line){  #seqs in different otus are separated by tab
			my @seqs=split(/,/,$otu);  #seqs in same otus are separated by ","
			$Fp+=@seqs*(@seqs-1)/2;  #how many pairs in this otu
			$Tn+=@seqs*($seqcount-@seqs)/2; #how many pairs between this and another otu
			foreach my $seq(@seqs){
				$otuassign{$seq}=$otuid; #seqs in the same otu have same $otuid
			}
		$otuid++;
		}
	}
}
close(LIST);
print "List file reading done!\n";

print "Reading dist file...\n";
open(DIST,$dist);
while(<DIST>){
	chomp;
	my @line=split(/\s/,$_);
	my $otua=$otuassign{$line[0]};
	my $otub=$otuassign{$line[1]};
	if($otua && $otub && $line[2]<=0.03){  #two seqs have similarity >=97%, should delete corresponding Fp and Tn
		if($otua eq $otub){
			$Fp--;
			$Tp++;
		}else{
			$Tn--;
			$Fn++;
		}
	}
} 

close(DIST);
print "Dist file reading done!\n";

$MCC=($Tp*$Tn-$Fp*$Fn)/sqrt(($Tp+$Fp)*($Tp+$Fn)*($Tn+$Fp)*($Tn+$Fn));

open(OUT,">$out");
print OUT "Tp\tTn\tFp\tFn\tMCC\n";
print OUT "$Tp\t$Tn\t$Fp\t$Fn\t$MCC";
print "Tp\tTn\tFp\tFn\tMCC\n";
print "$Tp\t$Tn\t$Fp\t$Fn\t$MCC\n";
