#!/usr/bin/perl -w

# For Qiime's pick_otu_through_otu_table.py script, When the py script prodeces rdp-otu-table file
# The rm_sub.pl can modify otu-table file to a more understandable format, and also summarizes hierarchy information

# the function is like ClasifyRichness.pl or PhyloRichness.pl but this script take more things into consideration
# and maybe better 

# Wang Guoyang, 2012-02-24

use Getopt::Long;

my $otutable = undef;

&GetOptions(
			'otutable|ot=s' => \$otutable,
			);

open TABLE,$otutable or die ("USAGE: rm_sub.pl --otutable/ot otutable.txt\n");

my @subs = qw/Acidimicrobidae Actinobacteridae Coriobacteridae Rubrobacteridae Sphaerobacteridae Acidimicrobineae Actinomycineae Actinopolysporineae Catenulisporineae Coriobacterineae Corynebacterineae Cystobacterineae Frankineae Glycomycineae Kineosporiineae Micrococcineae Micromonosporineae Nannocystineae Propionibacterineae Pseudonocardineae Rubrobacterineae Sorangiineae Sorangineae Sphaerobacterineae Streptomycineae Streptosporangineae/;

# my $len = @subs;

# print "totaly $len elements in array.\n"; # When the rdp-taxonomy updates, the former @subs line should update

my @notes = qw/TM7 SR1 OP10/; # These bacteria is special, they are directly hierarchied as phylum-genus

my %otuinfo;

my (%phylum,%class,%order,%family,%genus);

my $rich_sum = 0;

while (<TABLE>){
chomp;
if(/\#/){
next;
}
if(/^(\d+)\t(\d+)\tRoot\;([\S\s]+)/){
my $uid = $1;
my $urich = $2;
my $uphylo = $3;
$uphylo =~ s/"//g;
$uphylo =~ s/;/\t/g;
foreach my $sub(@subs){
		
		print STDOUT $sub." removed with otu $uid.\n" if ($uphylo =~ s/$sub\t//);
		}
print $1." removed with $uid\n" if($uphylo =~ s/(\S+ \d+)\t//);# eg. Carnobacteriaceae 2 removed with 5
my @up = split "\t",$uphylo;
if (@up<6){
		my $tem_end = $up[$#up];
		my $len = @up;
		$len = 6-$len;
		for(1..$len){
		push @up,"unclassified_".$tem_end;
		}

}

$uphylo = join("\t",@up);

foreach my $spe(@notes){
$uphylo = "$spe	$spe	$spe	$spe	$spe	$spe" if ($uphylo =~ /$spe/);
}

my @phy = split("\t",$uphylo);
$phylum{$phy[1]}="phylum".$phy[1];# add prefix to make every value unique, or the sum up result maybe strange in the next part
$class{$phy[2]}="class".$phy[2];
$order{$phy[3]}="order".$phy[3];
$family{$phy[4]}="family".$phy[4];
$genus{$phy[5]}="genus".$phy[5];


$otuinfo{$uid}->{Rich}=$urich; # No need for 2-d array or so-called hash-array, just use it the easiest way
$otuinfo{$uid}->{Phylo}=$uphylo;
$rich_sum+=$urich;

$otuinfo{$uid}->{Phylum} = $phy[1];
$otuinfo{$uid}->{Class} = $phy[2];
$otuinfo{$uid}->{Order} = $phy[3];
$otuinfo{$uid}->{Family} = $phy[4];
$otuinfo{$uid}->{Genus} = $phy[5];
}
}

close TABLE;

my $out = "rmsub_".$otutable;
open OUT,">$out";
select OUT;

## Sum up taxonomy-richness
foreach my $key(sort keys %phylum){
	foreach my $id(sort keys %otuinfo){
		$phylum{$key}->{Rich} += $otuinfo{$id}->{Rich} if ($key eq $otuinfo{$id}->{Phylum});
		}
		}
		
foreach my $key(sort keys %class){
	foreach my $id(sort keys %otuinfo){
		$class{$key}->{Rich} += $otuinfo{$id}->{Rich} if ($key eq $otuinfo{$id}->{Class});
		}
		}
foreach my $key(sort keys %order){
	foreach my $id(sort keys %otuinfo){
		$order{$key}->{Rich} += $otuinfo{$id}->{Rich} if ($key eq $otuinfo{$id}->{Order});
		}
		}
foreach my $key(sort keys %family){
	foreach my $id(sort keys %otuinfo){
		$family{$key}->{Rich} += $otuinfo{$id}->{Rich} if ($key eq $otuinfo{$id}->{Family});
		}
		}
foreach my $key(sort keys %genus){
	foreach my $id(sort keys %otuinfo){
		$genus{$key}->{Rich} += $otuinfo{$id}->{Rich} if ($key eq $otuinfo{$id}->{Genus});
		}
		}


print "Phylum\tNumber\tPercent\(\%\)\n";
my $check = 0;
foreach my $key(sort keys %phylum){
	$check += $phylum{$key}->{Rich};
	my $per = 100*($phylum{$key}->{Rich}/$rich_sum);
	printf "%s\t%d\t%0.4f\n",$key,$phylum{$key}->{Rich},$per;
	}
warn("phylum level summary not match the total number.with check is $check and summary is $rich_sum\n") if ($check ne $rich_sum);

print "\n\n";


print "Class\tNumber\tPercent\(\%\)\n";
$check = 0;
foreach my $key(sort keys %class){
	$check += $class{$key}->{Rich};
	my $per = 100*($class{$key}->{Rich}/$rich_sum);
	printf "%s\t%d\t%0.4f\n",$key,$class{$key}->{Rich},$per;
	}
warn("class level summary not match the total nuber.with check is $check and summary is $rich_sum\n") if ($check ne $rich_sum);

print "\n\n";

print "Order\tNumber\tPercent\(\%\)\n";
$check = 0;
foreach my $key(sort keys %order){
	$check += $order{$key}->{Rich};
	my $per = 100*($order{$key}->{Rich}/$rich_sum);
	printf "%s\t%d\t%0.4f\n",$key,$order{$key}->{Rich},$per;
	}
warn("order level summary not match the total nuber\n") if ($check ne $rich_sum);

print "\n\n";


print "Family\tNumber\tPercent\(\%\)\n";
$check = 0;
foreach my $key(sort keys %family){
	$check += $family{$key}->{Rich};
	my $per = 100*($family{$key}->{Rich}/$rich_sum);
	printf "%s\t%d\t%0.4f\n",$key,$family{$key}->{Rich},$per;
	}
warn("family level summary not match the total nuber\n") if ($check ne $rich_sum);

print "\n\n";


print "Genus\tNumber\tPercent\(\%\)\n";
$check = 0;
foreach my $key(sort keys %genus){
	$check += $genus{$key}->{Rich};
	my $per = 100*($genus{$key}->{Rich}/$rich_sum);
	printf "%s\t%d\t%0.4f\n",$key,$genus{$key}->{Rich},$per;
	}
warn("genus level summary not match the total nuber\n") if ($check ne $rich_sum);

print "\n\n\n";

print "\#otu ID\tNumber\tPercent\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\n";

foreach my $key(sort {$a<=>$b} keys %otuinfo){
my $percent = 100*($otuinfo{$key}->{Rich})/$rich_sum;
printf "%d\t%d\t%0.4f\t%s\n",$key,$otuinfo{$key}->{Rich},$percent,$otuinfo{$key}->{Phylo};
}





close OUT;
