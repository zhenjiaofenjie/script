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
# 程序第14行至20行，如果使用智能匹配或者不把$1,$2,$3在一开始就赋给其他变量，而是采用如下形式：
# $1 =~s/\-/\./g;
# $2 =~s/\-/\./g;
# print $1.$2.$3."\n"
# 这样是得不到结果的，或者只把$1,$3赋给一个其他变量，而后面使用 print $head.$2.$tail."\n";这样也是得不到结果的
# 智能匹配也不能直接做绑定后的替换，这可能是perl在匹配时的一些内部机制，而这种机制我现在并不知道
#=============================
