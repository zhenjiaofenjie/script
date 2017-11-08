#!/usr//bin/perl
use warnings;

my $usage="USAGE MergeDis.pl file1 file2\n";
$dis1 = $ARGV[0] or die($usage);
$dis2 = $ARGV[1] or die($usage);
open DIS1,$dis1;
open DIS2,$dis2;
my %dis1=();
my %dis2=();
my @array=();
while (<DIS1>){
if($.>1){
		chomp;
		@array = split;
		$dis1{$array[0]}=$array[1];
	}
}
my $lines1 = $.-1;
print "$lines1 UniV3ID in $dis1\n";
close DIS1;
while (<DIS2>){
if($.>1){
		chomp;
		@array = split;
		$dis2{$array[0]}=$array[1];
	}
}
my $lines2 = $.-1;
print "$lines2 UniV3ID in $dis2\n";
close DIS2;
foreach my $key(sort keys %dis1){
	if($dis2{$key}){
		$dis2{$key}+=$dis1{$key};
	}
	else {$dis2{$key}=$dis1{$key};}
}
open DIS,">merged.dis";
select DIS;
print "UniV3	×ÜÐòÁÐÊý\n";
foreach my $key(sort keys %dis2){
	print "$key	$dis2{$key}\n";
}
close DIS;
#========================================
# Written by Wang Guoyang
# 2011-04-24
# for meiging dis file after the extration work
# paul.paulo.king@gmail.com
#========================================