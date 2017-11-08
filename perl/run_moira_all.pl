#!/usr/bin/perl
use strict;
use warnings;

my $USAGE="perl run_moira_all.pl samplefile\n";
my $samplefile=shift or die($USAGE);
my $projectname=$samplefile;
$projectname=~ s/.txt//;
my $groupfile=$projectname."_moira.groups";
my $projectr1=$projectname.".R1.fastq";
my $projectr2=$projectname.".R2.fastq";

open(GROUP,">".$groupfile);
open(FASTQOUT1,">".$projectr1);   #To combine all original seq to one file.
open(FASTQOUT2,">".$projectr2);
open(FILE,$samplefile);
my $seqcount=1;
while(<FILE>){
	chomp;
	my @line=split(/\t/);
	my $samplename=shift(@line);
	my $fastqr1=shift(@line);
	my $fastqr2=shift(@line);
    open(FASTQ1,$fastqr1);
    open(FASTQ2,$fastqr2);
    while(<FASTQ1>){
        <FASTQ2>;
        print FASTQOUT1 "\@$samplename\_$seqcount;barcodelabel=$samplename;\n";        # To relabel the header of seqs
        print FASTQOUT2 "\@$samplename\_$seqcount;barcodelabel=$samplename;\n";
        print GROUP "$samplename\_$seqcount;barcodelabel=$samplename;\t$samplename\n";
        for(my $i=0;$i<3;$i++){ #print the rest of original seqs
            my $fastqin=<FASTQ1>;
            print FASTQOUT1 $fastqin;
            $fastqin=<FASTQ2>;
            print FASTQOUT2 $fastqin;
        }
        $seqcount++;
    }
    close(FASTQ1);
    close(FASTQ2);
} 
close(FASTQOUT1);
close(FASTQOUT2);
close(GROUP);
system("moira.py -ffq $projectr1 -rfq $projectr2 --paired -p 24");
system("mv $projectname\.R1.qc.good.fasta $projectname\_moira.fasta");
system("mv $projectname\.R1.qc.good.names $projectname\_moira.names");
