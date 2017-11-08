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
	@arrayforrich = split;my $otuhere = shift(@arrayforrich);$sum=0;#ʹ��shift֮�������Ѿ��仯����һ��Ԫ���Ƴ���������Ԫ�ر���Ӧ��0Ԫ�ؿ�ʼ
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
#����һ������ļ�����һ���ܵ�fasta�ļ���ѡ������ļ��г��ֵ�unique
#��������Щunique���µ�fasta�ļ�
#������������������е����ݽ��кϲ�����������ķ������������˲��˵����������Ժϲ�ʵ����֮ǰ��õ�����������
#�������ķ���������ɲο������������ҽԺ����A00041��A00001���ݵĻ��
###############################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#�����г�����һ������ϵĴ���,��ֱ��ʹ�ô˳�����г�ȡ������Ҫ����fasta�ļ��е����һ��Unique���е������ֹ�����
#||||||||||||||||||||||||||||||||||||||||||||||
#������ʹ��before���ݵķ��������ڴ���֮ǰ������������Ϣ����֮ǰ�Ĺ�ϣ�����������������key������ֱ�value����
#�⣬���ǵ����һ��Ԫ��ʱ������������û��������Ԫ�ع������ݣ���������һ��ֻ�м��Ŀչ�ϣ���Ӷ���������Ĵ���
#��ʱ��ʹ�� ����==��ϣ ��Ӧ�ķ��������ÿ���ǰ�ݻ��Ǻ��ݣ��ṹ�򵥣��Ƚϰ�ȫ��
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
################################################
# Written by Wang Guoyang,MOLECO,SJTU
# Email:paul.paulo.king@gmail.com
# contact 13761023529 if you run into any problem while using this script
# 2011-01-19
# 2011-03-15�޸ģ���������ɲ�����Ϊ0�з���ļ������
# 2011-05-03�޸ģ�������֮ǰ���ݷ��������Ĵ����޸�Ϊ��������ʽ
################################################
