#!/usr/bin/perl

#Used to allot Pacbio seqs to each sample.
#Written by WangJing 2012-05-17

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;


my($seqfile,$tagfile)=(undef,undef);

GetOptions('sq|s=s' => \$seqfile,
	  );

my $usage="USAGE: allot2sample.pl -s seqfile\n";
if(!$seqfile){
	die($usage);	
}

my %id;
                   
my $seq_in=Bio::SeqIO->new(-file=>$seqfile,
			   -format=>"fasta", 
                          );
                                                    
while(my $seq=$seq_in -> next_seq()){
	my $name=$seq->display_name();
	$name=~ s/_.*$//;
	push(@{$id{$name}},$seq);
}


foreach my $sample(keys(%id)){
	open(OUT,">".$sample.".fasta");
	while(${$id{$sample}}[0]){
		my $seq=shift(@{$id{$sample}});
		print OUT ">".$seq->display_id."\n".$seq->seq."\n";
	}
}
