#!/usr/bin/perl
use warnings;
use strict;

my $usage="pick_seed_from_uclust_uc.pl uclust.uc uniqmap.txt\n";
my $uc=shift or die($usage);
my $uniqmap=shift or die($usage);
my $seedout=$uc.".seed";
my $otuout=$uc.".otu.txt";

my %maphash;
open(MAP,$uniqmap);
while(<MAP>){
	chomp;
	my @line=split(/\t/);
	$maphash{shift @line}=join("\t",@line[1..@line-1]);
}
close(MAP);

my %seedhash;
my %otuhash;
open(UC,$uc);
while(<UC>){
	if(/^S/){   #new seed
		chomp;
		my @line=split(/\t/);
		$line[8]=~ /size=(\d*)/;
		$seedhash{$line[8]}+=$1;
		$otuhash{$line[8]}=$maphash{$line[8]};
	}elsif(/^H/){   #mapped to exist seed
		chomp;
		my @line=split(/\t/);
		$line[8]=~ /size=(\d*)/;
		$seedhash{$line[9]}+=$1;
		$otuhash{$line[9]}.="\t".$maphash{$line[8]};
	}
}
close(UC);

open(SEEDOUT,">".$seedout);
open(OTUOUT,">".$otuout);
my $otucount=1;
foreach my $key(keys %seedhash){
	print SEEDOUT "$key\t$seedhash{$key}\n";
	print OTUOUT "denovootu".$otucount."\t".$otuhash{$key}."\n";
	$otucount++;
}
close(SEEDOUT);
close(OTUOUT);
