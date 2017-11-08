#!/usr/bin/perl
use strict;
use warnings;

my $usage="sens_spec.pl listfile distfile output\n";
my $list=shift or die($usage);
my $dist=shift or die($usage);
my $out=shift or die($usage);

print "Reading dist file...\n";
open(DIST,$dist);
my %gooddist;
while(<DIST>){
    chomp;
    my @line=split(/\s/,$_);
    if($line[2]<=0.03){  #two seqs have similarity >=97%
        ${$gooddist{$line[0]}}{$line[1]}=1;
    }
}
close(DIST);
print "Dist file reading done!\n";

print "Reading list file...\n";
open(LIST,$list);
my %otumate;
my @seqnames;
while(<LIST>){
    if(/^0.03/){
        chomp;
        my @line=split(/\t/,$_);
        splice(@line,0,2);
        foreach my $otu(@line){  #seqs in different otus are separated by tab
            my @seqs=split(/,/,$otu);  #seqs in same otus are separated by ","
            push(@seqnames,@seqs);
            for(my $i=0;$i<@seqs-1;$i++){  #record the seqs in the same otu
                for(my $j=$i+1;$j<@seqs;$j++){
                    ${$otumate{$seqs[$i]}}{$seqs[$j]}=1;
                }
            }
        }
    }
}
close(LIST);
print "List file reading done!\n";

my $count;
my($Tp,$Tn,$Fp,$Fn,$MCC)=(0,0,0,0,0);
for(my $i=0;$i<@seqnames-1;$i++){  #for each pair of seqs
    for(my $j=$i+1;$j<@seqnames;$j++){
        $count++;
        if(${$gooddist{$seqnames[$i]}}{$seqnames[$j]} || ${$gooddist{$seqnames[$j]}}{$seqnames[$i]}){
            if(${$otumate{$seqnames[$i]}}{$seqnames[$j]} || ${$otumate{$seqnames[$j]}}{$seqnames[$i]}){
                $Tp++;
            }else{
                $Fn++;
            }
        }else{
            if(${$otumate{$seqnames[$i]}}{$seqnames[$j]} || ${$otumate{$seqnames[$j]}}{$seqnames[$i]}){
                $Fp++;
            }else{
                $Tn++;
            }
        }
    }
}
$MCC=($Tp*$Tn-$Fp*$Fn)/sqrt(($Tp+$Fp)*($Tp+$Fn)*($Tn+$Fp)*($Tn+$Fn));

open(OUT,">$out");
print OUT "Tp\tTn\tFp\tFn\tMCC\n";
print OUT "$Tp\t$Tn\t$Fp\t$Fn\t$MCC";
