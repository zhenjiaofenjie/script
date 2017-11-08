#!/usr/bin/perl
use strict;
use warnings;

my $USAGE="perl run_sickle+BH+pandaseq.pl samplefile\n";
my $samplefile=shift or die($USAGE);
my $projectname=$samplefile;
$projectname=~ s/\.txt//;

open(FILE,$samplefile);
while(<FILE>){
	chomp;
	my @line=split(/\t/);
	my $samplename=shift(@line);
	my $fastqr1=shift(@line);
	my $fastqr2=shift(@line);
	system("~/sickle-master/sickle pe -f $fastqr1 -r $fastqr2 -o $samplename\_R1_trimed.fastq -p $samplename\_R2_trimed.fastq -n -t sanger -s $samplename\_trimed_single.fastq");
 	system("~/SPAdes-3.5.0-Linux/bin/spades.py --careful --only-error-correction -1 $samplename\_R1_trimed.fastq -2 $samplename\_R2_trimed.fastq -o . --disable-gzip-output -t 12");
	system("sed '/^\@M03307/s/\$/ 1:N:0:/' corrected/$samplename\_R1_trimed.00.0_0.cor.fastq > $samplename\_R1_trimed_cor.fastq");
	system("sed '/^\@M03307/s/\$/ 2:N:0:/' corrected/$samplename\_R2_trimed.00.0_0.cor.fastq > $samplename\_R2_trimed_cor.fastq");
	system("pandaseq -B -f $samplename\_R1_trimed_cor.fastq -r $samplename\_R2_trimed_cor.fastq -o 50 -l 400 -w $samplename\_pandaseq.fasta -g $samplename\_pandaseq.log");
	open(SEDIN,"$samplename\_pandaseq.fasta");
	open(SEDOUT,">$samplename\_pandaseq_sed.fasta");
	my $count=1;
	while(<SEDIN>){
		chomp;
		if(/>/){	
			$_=">$samplename\_$count;barcodelabel=$samplename;";
			$count++;
		}
		print SEDOUT $_."\n";
	}
	close(SEDIN);
	close(SEDOUT);
}
