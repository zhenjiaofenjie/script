#!/usr/bin/perl -w

# ==========================================================================================================================================================
# This script is for remove RDP sub-classified hierarchies
# the @rdp_subs array contains sub-classifies of RDP hierarchy 2012-09-13
# getting from http://rdp.cme.msu.edu/hierarchy/hierarchy_browser.jsp?qvector=255&depth=10&openNode=0&seqid=&currentRoot=0&searchStr=&endDataValue=&showOpt=
# ==========================================================================================================================================================

# ===================================================================================
# The second application of the script is filter hierarchies under a given threashold
# when some hierarchy result in the phylo file is lower than the threashold
# the classified result is substituted by an "unclassified_NAME" item
# in which NAME is the upper level of the result
# ===================================================================================

# =======================
# Written by Wang Guoyang
# 2012-09-13 
# =======================

# UNFINISHED YET 

# Current implements sub-classifications removal and adding unclassified by recognizing
# threshold
# but the special case for TM7 liked unsolved, just don't forget checking the final result

use YAML;
use Getopt::Long;

my ($file, $threshold) = (undef, 0.8);

GetOptions(
			"f=s" => \$file,
			"th=f" => \$threshold,
			);


open PHYLO, $file or die("USAGE: RemoveSubtype.pl rdpfile.txt\n");

my @rdp_subs = qw(Acidimicrobidae Acidimicrobineae Actinobacteridae Actinomycineae Actinopolysporineae Catenulisporineae Corynebacterineae Frankineae Glycomycineae Jiangellineae Kineosporiineae Micrococcineae Micromonosporineae Propionibacterineae Pseudonocardineae Streptomycineae Streptosporangineae Coriobacteridae Coriobacterineae Nitriliruptoridae Rubrobacteridae Rubrobacterineae Sphaerobacteridae Sphaerobacterineae Cystobacterineae Nannocystineae Sorangiineae);

my %rdp_subs;

map {$rdp_subs{$_} = 1} @rdp_subs;

# print "Phylum\tClass\tOrder\tFamily\tGenus\n";

foreach (<PHYLO>){
	chomp;
	remove_subs_first();
}



sub form_line_hierarchy{
	
	# ============================================================
	# This function can be used for generating a line-based
	# tree, which element has 'parent', 'child', and 'similarity',
	# using a data structure of linked array to form a linked hash
	# But when the child has an identity name with the parent
	# such as Actinobacteria, the class name is just definitely
	# same as the phylum name, the function just cannot recognize..
	#
	# Anyway, I'll try to use a Module for this job and find the 
	# difference...
	#
	# Or, just using the classic Phylum..Genus method, it's ugly,
	# but it works..fine..also die sometime while sub-classificat-
	# -ions appear, on this case, removing all subs first will be
	# fine..I think, Ok, let's try
	
	# ============================================================
	
	my @line = split "\t",$_;
	my $name_no = ($#line - 1)/2 - 1;
	my $root = 'Bacteria';
	$line[0]->{Parent} = $root;
	$line[0]->{Child} = $line[2];
	$line[0]->{Similarity} = $line[1];
	map {$line[2*$_]->{Similarity} = $line[2*$_+1];
	$line[2*$_]->{Child} = $line[2*$_ + 2];
	$line[2*$_]->{Parent} = $line[2*$_ - 2];}
	1..$name_no;
	$line[$#line - 1]->{Parent} = $line[$#line - 3];
	$line[$#line - 1]->{Child} = 'None';
	$line[$#line - 1]->{Similarity} = $line[$#line];
	$line[0];
}


sub remove_subs_first{

	my $line = $_;
	map {$line =~ s/$_\t\d+\%\t//g} keys %rdp_subs;
	
#	print $line."\n";  ## for double check results with the original phylo line
	
	
	my @line = split "\t",$line;		
#	warn("some thing wrong:\n@line\n") if @line != 12;

#   ======================================================
#   lots of lines with subs will not having 12 elements...
#   so this cannot be a warning criteria
#   ======================================================
	
	my (@names,@similarities);
	
	my $name_no = ($#line - 1)/2;
	
	for (0..$name_no){		
		my $value  = $line[2*$_ + 1];
		$value =~ s/\%//;
		next if $value < ($threshold*100);
		$names[$_] = $line[2*$_];
		$similarities[$_] = $line[2*$_ + 1];
	}
	
	
	if (@names == 0){
		@names = qw(unclassified_bacteria unclassified_bacteria unclassified_bacteria unclassified_bacteria unclassified_bacteria);
	}
	
	my @for_names = @names;
	
	map {$names[$_] = "unclassified_".$for_names[$#for_names]} $#for_names+1..4 if $#for_names < 4;
	
	for (0..$#names){
#		my $heihei = defined($similarities[$_]) ? $similarities[$_] : 'heihei';
#		print $names[$_]."\t".$heihei."\t";
		
		print  $names[$_]."\t";
	}
	
	print "\n";
	
	# ============================================================
	# It's a pity that the rdp result is so strange,
	# something like TM7 have a instant genus after phylum
	# and after some subs removed, you cannot recognize some level
	# as upper level or lower level..
	# God help me..
	# ============================================================
}

sub correct_phylo_table{
	1; # TODO job
}
