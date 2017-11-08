#!/usr/bin/perl
use strict;
use warnings;
#To change the "size" annotation of filter.fasta based on the real count in filtered.count_table output by mothur.moira.step1.txt and mothur.moira.step2.txt.

my $usage="USAGE: resize_moira_filtered.pl filtered.fasta filtered.count_table\n";
my $fasta=shift or die($usage);
my $counttable=shift or die($usage);
my $fastaout=$fasta;
$fastaout=~ s/fasta/resize.fasta/;
my $fastanogap=$fastaout;
$fastanogap=~ s/fasta/nogap.fasta/;
my $countout=$counttable;
$countout=~ s/count_table/resize.count_table/;

open(COUNTTABLE,$counttable);
open(FASTA,$fasta);
open(FASTAOUT,">".$fastaout);
open(FASTANOGAP,">".$fastanogap);
open(COUNTOUT,">".$countout);

my %counthash;
my $_=<COUNTTABLE>;
print COUNTOUT;
while(<COUNTTABLE>){
    my @line=split(/\t/,$_);
    $counthash{$line[0]}=$line[1];
    if($_=~ /size=/){
        $_=~ s/size=\d*/size=$line[1]/;
    }else{
        $_=~ s/\t/size=$line[1];\t/;
    }
    print COUNTOUT $_;
}
close(COUNTTABLE);
close(COUNTOUT);

while(<FASTA>){
    if (/^>/){
        my $seqname=$_;
        chomp($seqname);
        $seqname=substr($seqname,1);
        if($_=~ /size=/){
            $_=~ s/size=\d*/size=$counthash{$seqname}/;
        }else{
            $_=~ s/$/size=$counthash{$seqname};/;
        }
        print FASTAOUT $_;
        print FASTANOGAP $_;
    }else{
        print FASTAOUT $_;
        $_=~ s/-//g;
        print FASTANOGAP $_;
    }
}
close(FASTA);
close(FASTAOUT);
close(FASTANOGAP);
