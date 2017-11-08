#!/usr/bin/perl
use strict;
use warnings;

my $usage="perl run_threshold.pl summed_threshold.txt\n";
my $file=shift or die($usage);
open(FILE,$file);
while(<FILE>){
	chomp;
	my @line=split(/\t/,$_);
	my $percent=$line[0];
	my $count=$line[1];
	system("mkdir unique_mc$count\_$percent/");
	system("~/usearch8 -sortbysize mock_usearch_pcrstriped_uniques.fasta -fastaout unique_mc$count\_$percent/mock_usearch_pcrstriped_uniques_sorted.fasta -minsize $count");
	system("~/usearch8 -cluster_otus unique_mc$count\_$percent/mock_usearch_pcrstriped_uniques_sorted.fasta -otus unique_mc$count\_$percent/mock_usearch_pcrstriped_otus.fasta");
	system("~/usearch8 -uchime_ref unique_mc$count\_$percent/mock_usearch_pcrstriped_otus.fasta -db ~/rdp_gold.fa -uchimeout unique_mc$count\_$percent/mock_usearch_pcrstriped_results.uchime -strand plus -chimeras unique_mc$count\_$percent/mock_usearch_pcrstriped_otus_chimeras.fasta");
	system("filter_fasta.py -f unique_mc$count\_$percent/mock_usearch_pcrstriped_otus.fasta -o unique_mc$count\_$percent/mock_usearch_pcrstriped_otus_nonchimeras.fasta -a unique_mc$count\_$percent/mock_usearch_pcrstriped_otus_chimeras.fasta -n");
	system("fasta_number.py unique_mc$count\_$percent/mock_usearch_pcrstriped_otus_nonchimeras.fasta OTU > unique_mc$count\_$percent/mock_usearch_pcrstriped_otus_numbered.fasta");
	system("~/usearch8 -usearch_global mock_usearch_pcrstriped.fasta -db unique_mc$count\_$percent/mock_usearch_pcrstriped_otus_numbered.fasta -strand plus -id 0.97 -uc unique_mc$count\_$percent/mock_usearch_pcrstriped_readmap.uc");
	system("python ~/.perl_program/uc2otutab.py unique_mc$count\_$percent/mock_usearch_pcrstriped_readmap.uc > unique_mc$count\_$percent/mock_usearch_pcrstriped_otu_table.txt");
	system("biom convert -i unique_mc$count\_$percent/mock_usearch_pcrstriped_otu_table.txt -o unique_mc$count\_$percent/mock_usearch_pcrstriped_otu_table.biom --table-type 'OTU table' --to-json");
	system("biom summarize-table -i unique_mc$count\_$percent/mock_usearch_pcrstriped_otu_table.biom -o unique_mc$count\_$percent/mock_usearch_pcrstriped_otu_table_summary.txt");
	system("pick_de_novo_otus.py -i unique_mc$count\_$percent/mock_usearch_pcrstriped_uniques_sorted.fasta -o unique_mc$count\_$percent/denovo -af");
	system("pick_closed_reference_otus.py -i mock_usearch_pcrstriped.fasta -r unique_mc$count\_$percent/denovo/rep_set/mock_usearch_pcrstriped_uniques_sorted_rep_set.fasta -o unique_mc$count\_$percent/denovo/closed/ -af");
	system("parallel_identify_chimeric_seqs.py -i unique_mc$count\_$percent/denovo/pynast_aligned_seqs/mock_usearch_pcrstriped_uniques_sorted_rep_set_aligned.fasta -o unique_mc$count\_$percent/denovo/chimeric.txt -a ~/qiime_software/core_set_aligned.fasta.imputed");
	system("filter_fasta.py -f unique_mc$count\_$percent/denovo/rep_set/mock_usearch_pcrstriped_uniques_sorted_rep_set.fasta -o unique_mc$count\_$percent/denovo/rep_set_nochimeric.fasta -s unique_mc$count\_$percent/denovo/chimeric.txt -n");
	system("~/usearch8 -uchime_ref unique_mc$count\_$percent/denovo/rep_set/mock_usearch_pcrstriped_uniques_sorted_rep_set.fasta -db ~/rdp_gold.fa -uchimeout unique_mc$count\_$percent/denovo/rep_set.uchime -strand plus -chimeras unique_mc$count\_$percent/denovo/rep_set_uchime_chimeras.fasta");
	system("filter_fasta.py -f unique_mc$count\_$percent/denovo/rep_set/mock_usearch_pcrstriped_uniques_sorted_rep_set.fasta -a unique_mc$count\_$percent/denovo/rep_set_uchime_chimeras.fasta -n -o unique_mc$count\_$percent/denovo/rep_set_uchime_nochimeric.fasta");
	system("pick_closed_reference_otus.py -i mock_usearch_pcrstriped.fasta -r unique_mc$count\_$percent/denovo/rep_set_nochimeric.fasta -o unique_mc$count\_$percent/denovo/nochimeric_closed/ -af");
	system("pick_closed_reference_otus.py -i mock_usearch_pcrstriped.fasta -r unique_mc$count\_$percent/denovo/rep_set_uchime_nochimeric.fasta -o unique_mc$count\_$percent/denovo/uchime_nochimeric_closed/ -af");
}
