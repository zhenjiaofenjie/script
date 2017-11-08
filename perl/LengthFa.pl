#!/usr/bin/perl -w

use Bio::SeqIO;
use Chart::Gnuplot;
use Getopt::Long;

my ($file,$step) = (undef,undef);

GetOptions ("fasta|fa=s" => \$file, 
                        "step|s=i"   => \$step, 
); 

my $usage="USAGE: LengthFa.pl --fasta/fa seqfile.fasta --step/s step_length\n";
open FA,$file or die $usage;
my $len = $file;
$len =~ s/\.[^.]+$/len.txt/;
my $jpg = $file;
$jpg =~ s/\.[^.]+$/len.jpg/;
open LEN,">$len";
select LEN;

my $fasbio = Bio::SeqIO->new(-file => $file,-format => 'FASTA');

while (my $lin = $fasbio->next_seq){
	print $lin->length."\n";
}

close LEN;
close FA;

my $line = 1;

select STDOUT;
open FILE,$len or die $usage;

my %data;
$line--; # Get value from this column, bu
while(<FILE>){
chomp;
my @arr = split;
$data{$arr[$line]}++;
}

my @values = sort {$a <=> $b} keys %data;

my $lowest = $values[0]/$step;
$lowest = int($lowest);

my $highest = $values[$#values]/$step;
$highest = int($highest);

$highest++;


my %grouped_data;

foreach my $val(@values){
my $ste = int($val/$step)*$step;
$grouped_data{$ste}+=$data{$val};
}

@values = ();

for my $ite($lowest..$highest){
$ite*=$step;
push @values,$ite;
}


my @data_values;

foreach my $va(@values){
if($grouped_data{$va}){
push @data_values,$grouped_data{$va};}
else{push @data_values,0;}
}


# Data
    my @x = @values;
    my @y = @data_values;
    my @data_sorted = @y;
@data_sorted = sort {$a <=> $b} @data_sorted;
print "@data_sorted"."\n";
    
    # Create chart object and specify the properties of the chart
    my $chart = Chart::Gnuplot->new(
        output => $jpg,
        title  => "Sequence length Distribution",
        xlabel => "Sequence length",
        ylabel => "Sequence Number",
	xrange => [$values[0],$values[$#values]],
	yrange => [0,1.25*$data_sorted[$#data_sorted]],
    );
    
    # Create dataset object and specify the properties of the dataset
    my $dataSet = Chart::Gnuplot::DataSet->new(
        xdata => \@x,
        ydata => \@y,
  #      title => "Plotting a line from Perl arrays",
        style => "boxes",
	fill =>1,
    );
    
    # Plot the data set on the chart
    $chart->plot2d($dataSet);
