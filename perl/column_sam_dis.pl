#!/usr/bin/perl
use strict;
use warnings;
use Excel::Writer::XLSX;
use Getopt::Long;

my($file,$step,$out)=
(  undef,undef,undef);

Getopt::Long::GetOptions(
	'file=s'	=>\$file,
	'step=s'	=>\$step,
	'out=s'		=>\$out,
);

my $usage="USAGE: column_sam_dis.pl -file sam_dis.txt -step number -out outfilename";
if(!$file||!$step||!$out){
	die($usage);
}

my $workbook  = Excel::Writer::XLSX->new( $out );
my $worksheet = $workbook->add_worksheet();
my $bold      = $workbook->add_format( bold => 1 );

# Add the worksheet data that the charts will refer to.
my %hash;
open (FILE,$file);
my $max=0;
while(<FILE>){
	chomp;
	my @value=split(/\t/,$_);
	$hash{int($value[1]/$step)}++;
	if($max<int($value[1]/$step)){
		$max=int($value[1]/$step);
	}
}
close(FILE);

for(my $row=0;$row<=$max;$row++){
	my @value=($row*$step."-".($row+1)*$step,$hash{$row});
	$worksheet->write($row,0,\@value);
}

# Create a new chart object. In this case an embedded chart.
my $chart1 = $workbook->add_chart( type => 'column', embedded => 1 );

$chart1->add_series(
    categories => [ 'Sheet1', 0, $max, 0, 0 ],
    values     => [ 'Sheet1', 0, $max, 1, 1 ],
);

# Add a chart title and some axis labels.
$chart1->set_title ( name => 'Extract Distribution' );
$chart1->set_x_axis( name => 'Number of Seqs' );
$chart1->set_y_axis( name => 'Number of Samples' );

# Insert the chart into the worksheet (with an offset).
$worksheet->insert_chart( 'D2', $chart1);
