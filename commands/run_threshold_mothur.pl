#!/usr/bin/perl
use strict;
use warnings;

my $usage="perl run_threshold.pl summed_threshold.txt\n";
my $file=shift or die($usage);
system("convert_mothur_fasta_for_qiime.pl GD.contigs.good.pcr.groups GD.trim.contigs.good.pcr.fasta GD.trim.contigs.good.pcr.qiime.fasta");
open(FILE,$file);
while(<FILE>){
	chomp;
	my @line=split(/\t/,$_);
	my $percent=$line[0];
	my $count=$line[1];
	system("mkdir unique_mc$count\_$percent/");
	system("cp GD.trim.contigs.good.pcr.unique.good.filter.unique.fasta GD.trim.contigs.good.pcr.unique.good.filter.count_table mothur.command.threshold.txt unique_mc$count\_$percent/");
	chdir("unique_mc$count\_$percent/");
	system("sed -i '/split.abund/s/cutoff=\\w*/cutoff=$count/' mothur.command.threshold.txt");
	system("~/mothur/mothur mothur.command.threshold.txt");
	#system("sed -e '/^>/!s/-//g' -e '/^>/s/.*\\(Otu[0-9]*\\).*/>\\1/' GD.trim.contigs.good.pcr.unique.good.filter.unique.abund.precluster.pick.pick.an.unique_list.0.03.rep.fasta >GD.an.0.03.rep.sed.fasta");
	system("pick_closed_reference_otus.py -i ../GD.trim.contigs.good.pcr.qiime.fasta -r GD.an.0.03.rep.sed.fasta -o closed/ -af");
	chdir("..");
}
