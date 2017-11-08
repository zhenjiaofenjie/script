#!/usr/bin/perl -w
use Bio::SeqIO;
use YAML;

my @fasta = @ARGV or die("USAGE: mergeFa.pl fasta1 fasta2 ...\n");
&mergefa(@fasta);


sub mergefa{
	my %fas_all;
	foreach(@_){
		my $fas = Bio::SeqIO->new(-file => $_, -format => 'fasta');
		while(my $seq = $fas->next_seq){
			$fas_all{$seq->display_id}->{DESC} = $seq->desc;
			$fas_all{$seq->display_id}->{SEQ} = $seq->seq;
		}
	}
	foreach my $key(sort keys %fas_all){
		&print_fasta($key,$fas_all{$key}->{DESC},$fas_all{$key}->{SEQ});
		}
	#print Dump(%fas_all);
}

sub print_fasta{
	my $sample = shift;
	my $desc = shift;
	my $sequence_fas = shift;
	print ">".$sample." ".$desc."\n".$sequence_fas."\n";
}
