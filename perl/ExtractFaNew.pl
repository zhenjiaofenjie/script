#!/usr/bin/perl
use warnings;
use 5.010;
my %unique;
my $sum;

$USAGE = "\nExtractFaNew.pl abundance.txt original.fa";

$abundance = $ARGV[0];

if (!defined $abundance) {
    die($USAGE . "\nPlease specify the two arguments.\n");
}

my $naming = "";
if ($abundance =~ /(\S+)\.txt/){
	$naming = $1;
}
else {
	die("$USAGE\nThe first argument should be a file with suffix of \".txt\" ,means a abundance file.\n");
	}

$fasta = $ARGV[1];
my $fa_naming = "";
if ($fasta =~ /(\S+)\.fa/){
	$fa_naming = $1;
}
else {
	die("$USAGE\nThe second argument should be a file with suffix of \".fa\" ,means a fasta file.\n");
	}

my $new_dis = $naming."_Real.dis";
my $new_fa = $naming."_Real.fa";
my $new_abun = $naming."_Real.txt";

open Rich,$abundance;
open NewRich,">$new_abun";
select NewRich;
while (<Rich>){
if($.==1){print;}
else {
	chomp;
	@arrayforrich = split;my $otuhere = shift(@arrayforrich);$sum=0;#使用shift之后数组已经变化，第一个元素移出，对数字元素遍历应从0元素开始
	for (my $ite=0;$ite<=$#arrayforrich;$ite++){
		$sum += $arrayforrich[$ite];
		}
		$unique{$otuhere}=$sum;
#		print "$otuhere ${$unique{$otuhere}}[0]\n";
if ($sum){print;print"\n";}
	}
}
close Rich;
close NewRich;
open DIS,">$new_dis";
select DIS;
foreach my $key(sort keys %unique){
if($unique{$key}){
	print "$key	$unique{$key}\n";}
}
close DIS;
my %fullunique;
open Fasta,$fasta;
my $unique_here = "";
while (<Fasta>){
	chomp;
	if(/^\>(?<ID>\S+) (?<Forprint>[\w=]+)/){
		$unique_here = $+{ID};
		${$fullunique{$unique_here}}[0] =$+{Forprint};
	}
	else {${$fullunique{$unique_here}}[1].=$_;}
}
close Fasta;
open Health,">$new_fa";
select Health;
foreach my $key(sort keys %unique){
	if($unique{$key}){
		if(${$fullunique{$key}}[0]){
	print "\>$key ${$fullunique{$key}}[0] num=$unique{$key}\n";
print ${$fullunique{$key}}[1]."\n";}
	else {print "This Unique $key does not exist in the original fasta.\n"}
}
#	else {print "There are 0 Unique $key summed from the abundance file.\n";}
}
close Health;
###############################################
#给定一个丰度文件，从一个总的fasta文件中选出丰度文件中出现的unique
#并作出这些unique的新的fasta文件
#这可用于新数据与已有的数据进行合并，以做更多的分析，比如测得了病人的样本，可以合并实验室之前测得的正常人样本
#做后续的分析，具体可参考王敬翰的瑞金医院病例A00041与A00001数据的混合
###############################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#程序中出现了一个设计上的错误,若直接使用此程序进行抽取，则需要对总fasta文件中的最后一个Unique进行单独的手工操作
#||||||||||||||||||||||||||||||||||||||||||||||
#程序中使用before回溯的方法，将内存中之前读到的序列信息赋给之前的哈希键，这样解决了先有key，后读分别value的问
#题，但是到最后一个元素时，由于它后面没有其他的元素供它回溯，所以这是一个只有键的空哈希，从而造成其他的错误。
#这时，使用 数组==哈希 对应的方法，不用考虑前溯还是后溯，结构简单，比较安全。
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
################################################
# Written by Wang Guoyang,MOLECO,SJTU
# Email:paul.paulo.king@gmail.com
# contact 13761023529 if you run into any problem while using this script
# 2011-01-19
# 2011-03-15修改，添加了生成不含和为0行丰度文件的语句
# 2011-05-03修改，修正了之前回溯方法产生的错误，修改为两参数形式
################################################
