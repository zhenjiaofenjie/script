#!/usr/bin/perl -w
# Wang Guoyang 2012-02-16

if(!@ARGV){
	die("USAGE: merge_all_abund.pl file1 file2 ...\n");
}
use YAML;


my %hash2d;
my %key_dim2;


## READ HASH AND MERGE HASHES
foreach(@ARGV){
&read_2d_hash($_);
}


## PRINT MERGED HASHES
my @key2s = keys(%key_dim2);
print "sample";

foreach my $key2(@key2s){
print "\t".$key2;
}

print "\n";

foreach my $key(sort keys %hash2d){
print $key;
foreach my $key2(@key2s){
if(defined($hash2d{$key}{$key2})){
print "\t".$hash2d{$key}{$key2};
}
else {print "\t0";}
}
print "\n";
}



sub read_2d_hash{

open MAT,$_[0] or die $0;
my @keys;
my @values;

while (<MAT>){

chomp;
my $line = $_;

if($.==1){
$line =~ s/^(\S)*?\t//;
@keys = split("\t",$line);
foreach (@keys){
$key_dim2{$_} = 1;
}#end foreach
}#end if

else{
@values = split("\t",$line);
my $dir = shift @values;

warn "matrix not right." if($#keys!=$#values);

for my $inc(0..$#keys){
$hash2d{$dir}{$keys[$inc]} += $values[$inc];
}#end for

}#end else
}#end while
close MAT;
%hash2d;
}#end sub


