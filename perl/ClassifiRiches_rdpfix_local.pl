#!/usr/bin/perl
#===================================================
#将RDP Classifier的分类地位文件结合丰度文件，按分类地位重新统计丰度
#===================================================
#***************************************************
#Written by Wang Guoyang, 2011-01-05,MOLECO,SJTU
#contact 13761023529 if you run into any problem while using this script
# Email:paul.paulo.king@gmail.com
#***************************************************
#Distinctly modified by Jing Wang, 20140930
#wjingsjtu@gmail.com
#****************************************************
use warnings;
use strict;
my $phylofile = $ARGV[0];
my $richfile = $ARGV[1];
my $threshold=0.8;
if(@ARGV==3){
	$threshold=$ARGV[2];
}

my $USAGE = "USAGE:\nperl PhyloRichness_rdpfix_local.pl phylo_rdp_fixrank.txt abundance.txt <optional>threshold\n\nthreshold is the confidence threshold comes from RDP results, default 0.8\n";

if (!defined $phylofile) {
    die($USAGE . "\nPlease specify the two arguments.\n");
}

my $naming = "";
if ($richfile =~ /(\S+)\.txt/){
	$naming = $1;
}else {
	die("$USAGE\nThe Second argument should be a file with suffix of \".txt\" ,means a abundance file.\n");
}
my $phyloout = $naming."_phylo.txt";
my $otuout=$naming."_otuclassify.txt";
my (%kingdomname,%phylumname,%classname,%ordername,%familyname,%genusname);
open Phylofile,$phylofile;#此处为统计地位文件
open Richness,$richfile;#此处为丰度文件
foreach (<Phylofile>){
	chomp;
	$_=~ s/\"//g;
	if(/(?<UniqueID>[0-9a-zA-Z\#\.\_]+)\t[\-\+]*\t(?<Kingdom>[^\t]+)\tdomain\t(?<Kconfidence>[\d\.]+)\t(?<Phylum>[^\t]+)\tphylum\t(?<Pconfidence>[\d\.]+)\t(?<Class>[^\t]+)\tclass\t(?<Cconfidence>[\d\.]+)\t(?<Order>[^\t]+)\torder\t(?<Oconfidence>[\d\.]+)\t(?<Family>[^\t]+)\tfamily\t(?<Fconfidence>[\d\.]+)\t(?<Genus>[^\t]+)\tgenus\t(?<Gconfidence>[\d\.]+)/){
		my $kingdom=$+{Kingdom};
		my $phylum=$+{Phylum};
		my $class=$+{Class};
		my $order=$+{Order};
		my $family=$+{Family};
		my $genus=$+{Genus};
		if($+{Pconfidence}<$threshold){
		        $phylum=$kingdom;
			$phylum=~ s/Unclassified\s//g;
			$phylum="Unclassified ".$phylum;
		}
		if($+{Cconfidence}<$threshold){
			$class=$phylum;
			$class=~ s/Unclassified\s//g;
			$class="Unclassified ".$class;
		}
		if($+{Oconfidence}<$threshold){
			$order=$class;
                        $order=~ s/Unclassified\s//g;
			$order="Unclassified ".$order;
		}
		if($+{Fconfidence}<$threshold){
			$family=$order;
                        $family=~ s/Unclassified\s//g;
			$family="Unclassified ".$family;
		}
		if($+{Gconfidence}<$threshold){
			$genus=$family;
                        $genus=~ s/Unclassified\s//g;
			$genus="Unclassified ".$genus;
		}

		#print "$kingdom,$phylum,$class,$order,$family,$genus\n";
		#print $+{Gconfidence}."\n";
		
		$kingdomname{$+{UniqueID}}=$kingdom;
		$phylumname{$+{UniqueID}}=$phylum;
		$classname{$+{UniqueID}}=$class;
		$ordername{$+{UniqueID}}=$order;
		$familyname{$+{UniqueID}}=$family;
		$genusname{$+{UniqueID}}=$genus;
	}
}
my $linehead = "";
my (%kingdomrich,%phylumrich,%classrich,%orderrich,%familyrich,%genusrich);
open(OTUOUT,">".$otuout);
print OTUOUT "OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\n";
while (<Richness>){
	if(/^\#?OTU/){
		chomp;
		$linehead = $_;
	}elsif(/^\#/){
		next;
	}else {
		chomp;
		my @arrayforrich = split;
		my $otuid = shift(@arrayforrich);
		print OTUOUT $otuid."\t".$kingdomname{$otuid}."\t".$phylumname{$otuid}."\t".$classname{$otuid}."\t".$ordername{$otuid}."\t".$familyname{$otuid}."\t".$genusname{$otuid}."\n";
		for(my $i=0;$i<@arrayforrich;$i++){
			${$kingdomrich{$kingdomname{$otuid}}}[$i]+=$arrayforrich[$i];
			${$phylumrich{$phylumname{$otuid}}}[$i]+=$arrayforrich[$i];
			${$classrich{$classname{$otuid}}}[$i]+=$arrayforrich[$i];
			${$orderrich{$ordername{$otuid}}}[$i]+=$arrayforrich[$i];
			${$familyrich{$familyname{$otuid}}}[$i]+=$arrayforrich[$i];
			${$genusrich{$genusname{$otuid}}}[$i]+=$arrayforrich[$i];
		}
	}
}

open OUT,">$phyloout";
select OUT;
print "kingdom$linehead\n";
foreach my $key(sort keys %kingdomrich){
	print $key."\t".join("\t",@{$kingdomrich{$key}})."\n";
}
print "\nphylum$linehead\n";
foreach my $key(sort keys %phylumrich){
	print $key."\t".join("\t",@{$phylumrich{$key}})."\n";
}
print "\nclass$linehead\n";
foreach my $key(sort keys %classrich){
	print $key."\t".join("\t",@{$classrich{$key}})."\n";
}
print "\norder$linehead\n";
foreach my $key(sort keys %orderrich){
	print $key."\t".join("\t",@{$orderrich{$key}})."\n";
}
print "\nfamily$linehead\n";
foreach my $key(sort keys %familyrich){
	print $key."\t".join("\t",@{$familyrich{$key}})."\n";
}
print "\ngenus$linehead\n";
foreach my $key(sort keys %genusrich){
	print $key."\t".join("\t",@{$genusrich{$key}})."\n";
}
