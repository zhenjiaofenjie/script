#!/usr/bin/perl
#===================================================
#将RDP Classifier的分类地位文件结合丰度文件，按分类地位重新统计丰度
#使用时将此脚本放在上述两个文件所在目录下，修改脚本中对应的文件名
#生成的phylorich.txt即为最终的统计结果
#===================================================
#***************************************************
#Written by Wang Guoyang, 2011-01-05,MOLECO,SJTU
#contact 13761023529 if you run into any problem while using this script
# Email:paul.paulo.king@gmail.com
#***************************************************

#use warnings;
use 5.010;
my %phyloall;
my $phylofile = $ARGV[0];
my $richfile = $ARGV[1];

$USAGE = "USAGE:\nperl PhyloRichness.pl phylo.txt abundance.txt\n";

if (!defined $phylofile) {
    die($USAGE . "\nPlease specify the two arguments.\n");
}

my $naming = "";
if ($richfile =~ /(\S+)\.txt/){
	$naming = $1;
}
else {
	die("$USAGE\nThe Second argument should be a file with suffix of \".txt\" ,means a abundance file.\n");
	}
my $phyloout = $naming."_phylo.txt";
open Phylofile,$phylofile;#此处为统计地位文件
open Richness,$richfile;#此处为丰度文件
foreach (<Phylofile>){
chomp;
	if(/(?<UniqueID>[0-9a-zA-Z#]+);\+?\-?;Root;100\%;Bacteria;\d+\%;\"?\"?(?<Phylum>[a-zA-Z]+)\"?\"?;\d+\%;\"?\"?(?<Class>[a-zA-Z]+)\"?\"?;\d+\%;\"?\"?(?<Order>[a-zA-Z]+)\"?\"?;\d+\%;\"?\"?(?<Family>[a-zA-Z]+)\"?\"?;\d+\%;\"?\"?(?<Genus>[a-zA-Z]+)\"?\"?/){
#print "yes\n";
#print "$+{Phylum},$+{Class},$+{Order},$+{Family},$+{Genus}\n";
#print "$+{Species}\n";
#print "OK!\n";
push @{${$phyloall{$+{UniqueID}}}[0]},($+{Phylum},$+{Class},$+{Order},$+{Family},$+{Genus});}
}
$linehead = "";
$samplenum = 0;
my @arrayforrich=();
while (<Richness>){
if($.==1){chomp;$linehead = $_;$samplenum = s/\t/\t/g;}
else {chomp;@arrayforrich = split;my $otuhere = shift(@arrayforrich);push @{$phyloall{$otuhere}[1]},@arrayforrich;}
}
my (%phylum,%class,%order,%family,%genus);
open OUT,">$phyloout";
select OUT;
print "phylum$linehead\n";
foreach my $key(sort keys %phyloall){
#print "$key	${${$phyloall{$key}}[0]}[0] dododo\n";
push @{${$phylum{${${$phyloall{$key}}[0]}[0]}}[0]},$key;
if (${$phylum{${${$phyloall{$key}}[0]}[0]}}[1]){
for (my $item = 0;$item<$samplenum;$item++){
${${$phylum{${${$phyloall{$key}}[0]}[0]}}[1]}[$item]+=${${$phyloall{$key}}[1]}[$item];
}}
else {push @{${$phylum{${${$phyloall{$key}}[0]}[0]}}[1]},@{${$phyloall{$key}}[1]};}
}
foreach my $key(sort keys %phylum){
my $richhere = join "	",@{${$phylum{$key}}[1]};
print "$key	$richhere\n";
}
print "\nclass$linehead\n";
foreach my $key(sort keys %phyloall){
push @{${$class{${${$phyloall{$key}}[0]}[1]}}[0]},$key;
if (${$class{${${$phyloall{$key}}[0]}[1]}}[1]){
for (my $item = 0;$item<$samplenum;$item++){
${${$class{${${$phyloall{$key}}[0]}[1]}}[1]}[$item]+=${${$phyloall{$key}}[1]}[$item];
}}
else {push @{${$class{${${$phyloall{$key}}[0]}[1]}}[1]},@{${$phyloall{$key}}[1]};}
}
foreach my $key(sort keys %class){
my $richhere = join "	",@{${$class{$key}}[1]};
print "$key	$richhere\n";
}
print "\norder$linehead\n";
foreach my $key(sort keys %phyloall){
push @{${$order{${${$phyloall{$key}}[0]}[2]}}[0]},$key;
if (${$order{${${$phyloall{$key}}[0]}[2]}}[1]){
for (my $item = 0;$item<$samplenum;$item++){
${${$order{${${$phyloall{$key}}[0]}[2]}}[1]}[$item]+=${${$phyloall{$key}}[1]}[$item];
}}
else {push @{${$order{${${$phyloall{$key}}[0]}[2]}}[1]},@{${$phyloall{$key}}[1]};}
}
foreach my $key(sort keys %order){
my $richhere = join "	",@{${$order{$key}}[1]};
print "$key	$richhere\n";
}
print "\nfamily$linehead\n";
foreach my $key(sort keys %phyloall){
push @{${$family{${${$phyloall{$key}}[0]}[3]}}[0]},$key;
if (${$family{${${$phyloall{$key}}[0]}[3]}}[1]){
for (my $item = 0;$item<$samplenum;$item++){
${${$family{${${$phyloall{$key}}[0]}[3]}}[1]}[$item]+=${${$phyloall{$key}}[1]}[$item];
}}
else {push @{${$family{${${$phyloall{$key}}[0]}[3]}}[1]},@{${$phyloall{$key}}[1]};}
}
foreach my $key(sort keys %family){
my $richhere = join "	",@{${$family{$key}}[1]};
print "$key	$richhere\n";
}
print "\ngenus$linehead\n";
foreach my $key(sort keys %phyloall){
push @{${$genus{${${$phyloall{$key}}[0]}[4]}}[0]},$key;
if (${$genus{${${$phyloall{$key}}[0]}[4]}}[1]){
for (my $item = 0;$item<$samplenum;$item++){
${${$genus{${${$phyloall{$key}}[0]}[4]}}[1]}[$item]+=${${$phyloall{$key}}[1]}[$item];
}}
else {push @{${$genus{${${$phyloall{$key}}[0]}[4]}}[1]},@{${$phyloall{$key}}[1]};}
}
foreach my $key(sort keys %genus){
my $richhere = join "	",@{${$genus{$key}}[1]};
print "$key	$richhere\n";
}
