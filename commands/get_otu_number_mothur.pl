#!/usr/bin/perl
use strict;
use warnings;

#system("for i in unique*%/*.biom;do echo \$i;biom summarize-table -i \$i | head -n 3;done > usearch_otu.txt;");
system("for i in unique*%/closed/*.biom;do echo \$i;biom summarize-table -i \$i | head -n 3;done > mothur_otu.txt");
#system("for i in unique*%/denovo/nochimeric_closed/*.biom;do echo \$i;biom summarize-table -i \$i | head -n 3;done > denovo_nochimeric_otu.txt");
#system("for i in unique*%/denovo/uchime_nochimeric_closed/*.biom;do echo \$i;biom summarize-table -i \$i | head -n 3;done > denovo_uchime_nochimeric_otu.txt");
my @file=("mothur_otu.txt");
while(my $file=shift(@file)){
	open(FILE,$file);
	my $out=$file;
	$out=~s/\.txt/_summary.txt/;
	open(OUT,">".$out);
	print OUT "MinCount\t%discard\tNo.ofsample\tNo.ofOTU\tNo.ofseqs";
	while(<FILE>){
		if(/mc(\d*)_(\d*)\%/){
			print OUT "\n".$1."\t".$2;
		}else{
			chomp;
			my @line=split(/:\s/,$_);
			print OUT "\t".$line[1];
		}
	}
	close(FILE);
	close(OUT);
}
