#!/usr/bin/perl
use warnings;
use 5.010;
my $usage='USAGE: convert_pynast_to_arb.pl rep_seq_file';
my $with_file=$ARGV[0] or die($usage);
open FILE,$with_file;
open OUT,">$with_file.dot.fasta";
select OUT;
foreach (<FILE>){
chomp;
if (/^\>/){
	print $_."\n";
}
else {
if (/(\-+)([ATCGNRYKMSW][-ATCGNRYKMSW]+[ATCGNRYKMSW])(\-+)/){
my $head = $1;
my $mid = $2;
my $tail = $3;
$head =~ s/\-/\./g;
$tail =~ s/\-/\./g;
print $head.$mid.$tail."\n";
}
}
}
close FILE;
close OUT;
#=============================
# Written by Wang Guoyang
# 2011-04-01
# paul.paulo.king@gmail.com
# �����14����20�У����ʹ������ƥ����߲���$1,$2,$3��һ��ʼ�͸����������������ǲ���������ʽ��
# $1 =~s/\-/\./g;
# $2 =~s/\-/\./g;
# print $1.$2.$3."\n"
# �����ǵò�������ģ�����ֻ��$1,$3����һ������������������ʹ�� print $head.$2.$tail."\n";����Ҳ�ǵò��������
# ����ƥ��Ҳ����ֱ�����󶨺���滻���������perl��ƥ��ʱ��һЩ�ڲ����ƣ������ֻ��������ڲ���֪��
#=============================
