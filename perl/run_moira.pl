#!/usr/bin/perl
use strict;
use warnings;

my $USAGE="perl run_moira.pl samplefile\n";
my $samplefile=shift or die($USAGE);

open(FILE,$samplefile);
while(<FILE>){
	chomp;
	my @line=split(/\t/);
	my $samplename=shift(@line);
	my $fastqr1=shift(@line);
	my $fastqr2=shift(@line);
    system("moira.py -ffq $fastqr1 -rfq $fastqr2 --paired -p 24 -l $samplename\_");
    my $namefile=$fastqr1;
    $namefile=~ s/fastq/qc.good.names/;
    my $nameout=$namefile;
    $nameout=~ s/names/relabel.names/;
    my $groupfile=$namefile;
    $groupfile=~ s/names/groups/;
    my $fastafile=$fastqr1;
    $fastafile=~ s/fastq/qc.good.fasta/;
    my $fastaout=$fastafile;
    $fastaout=~ s/.fasta/.relabel.fasta/;
    open(NAMEFILE,$namefile);
    open(NAMEOUT,">".$nameout);
    open(GROUPFILE,">".$groupfile);
    open(FASTAIN,$fastafile);
    open(FASTAOUT,">".$fastaout); 
    while(<NAMEFILE>){
        chomp;
        my @line=split(/\t/,$_);
        my @seqnames=split(/,/,$line[1]);
        foreach my $seqname(@seqnames){
            print GROUPFILE "$seqname\t$samplename\n";
        }
        my $seqnum=@seqnames;
        my $fastaread=<FASTAIN>;
        chomp($fastaread);
        print FASTAOUT $fastaread.";barcodelabel=$samplename;size=$seqnum;\n";
        $fastaread=<FASTAIN>;
        print FASTAOUT $fastaread;
        print NAMEOUT $line[0].";barcodelabel=$samplename;size=$seqnum;\t".$line[1]."\n";
    }
    close(NAMEFILE);
    close(GROUPFILE);
    close(FASTAIN);
    close(FASTAOUT); 
} 
