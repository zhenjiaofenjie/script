#!/usr/bin/perl -w

# Wang Guoyang, 2012-08-29

if(!@ARGV){
	die("USAGE: fastq2fastaQual.pl fastqfile.fastq\n");
}
my $fastq = $ARGV[0];
open FASTQ,$ARGV[0] or die $!;

$fastq =~ s/(fastq|fq)//;

my $fasta = $fastq."fa";
my $qual = $fastq."qual";

open FASTA,">$fasta";
open QUAL,">$qual";

while (<FASTQ>){
	if($.%4==1){
		s/.*?-/>/;
		s/^@/>/;
		s/ .*//;
		print FASTA;
		print QUAL;
	}
	elsif($.%4==2){
		print FASTA;
	}
	elsif($.%4==0){
		my @q_c = qual_convert($_);
		print QUAL "@q_c\n";
	}
}

sub qual_convert{
	chomp;
	my @quals = split("");
	my @qual = map {ord($_)-33} @quals;
}
