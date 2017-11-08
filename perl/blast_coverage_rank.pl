#!/usr/bin/perl
# ===================================================
# 2011-12-05
# Show Blast results in the coverage-identity format
# ===================================================
# First time understands the 'Autovivification' concept
# Wang Guoyang
# ===================================================

use strict;
use Bio::SearchIO;

my $blastout = shift or die("USAGE: blast_coverage_rank.pl blastfile.out outfile\n");
my $out=shift or die("USAGE: blast_coverage_rank.pl blastfile.out outfile\n");
my ($report,$result,$hit,$rank);
my %res_cover;

open FOR, $blastout or die;
$report = Bio::SearchIO->new(-fh=>\*FOR, -format=>'blast');
while($result = $report->next_result){  # $result is a Bio::Search::Result::ResultI object
		while ($hit = $result->next_hit()) # # $hit is a Bio::Search::Hit::HitI object
		{while( my $hsp = $hit->next_hsp ) {
			my $coverage = ($hsp->length('hit'))/($result->query_length());
			unless(defined($res_cover{$result->query_name})){
		   $res_cover{$result->query_name}->{COVERAGE} = $coverage;
		   $res_cover{$result->query_name}->{IDENTITY} = $hsp->frac_identical('hit');
		   $res_cover{$result->query_name}->{DESCRIPTION} = $hit->description;
		   }
			elsif($res_cover{$result->query_name}<$coverage){
		   $res_cover{$result->query_name}->{COVERAGE} = $coverage;
		   $res_cover{$result->query_name}->{IDENTITY} = $hsp->frac_identical('hit');
		 $res_cover{$result->query_name}->{DESCRIPTION} = $hit->description;}
		}
	}
}

print STDOUT "blast file reading ok!\n";

open OUT,">$out";
select OUT;

print "OTUNAME	COVERAGE	IDENTITY	DESCRIPTION\n";
foreach my $key(sort keys %res_cover){
	print $key."	".$res_cover{$key}->{COVERAGE}."	".$res_cover{$key}->{IDENTITY}."	".$res_cover{$key}->{DESCRIPTION}."\n";
}

close OUT;