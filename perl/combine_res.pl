#!/usr/bin/perl
## Combine the Excel results produced by composition_excel.pl
## Input xlsx file numbers unlimited, all results gatheres
## Written by Wang Guoyang
## 2011-11-05
## paul.paulo.king@gmail.com

## ADVICE ****
## As the composition.pl exists some bugs, the xlsx used is STRONGLY RECOMMANDED for a check
## Especiely for the 'Incertae Sedis XIV', cause it's name contains whitespace

use Excel::Writer::XLSX;
use Spreadsheet::XLSX;
use YAML;
use Data::Dumper;
use strict;

#read_files();
#form_hashes();
#output_res();

## transverse xlsxs, each time open a workbook, form hashes for each separate file
## push hash references to an array
## transverse array, form a new hash, transverse this new hash, write to a new xlsx
my $usage="Combine the Excel results produced by composition_excel.pl\nUSAGE: combine_res.pl file1 file2 ...\n";
my @filelist = @ARGV or die($usage);
my $file_merged;
my %data;

foreach my $file(@filelist){
 my $excel = Spreadsheet::XLSX -> new ($file);
 my $fii = $file;
 $fii =~ s/_composition\.xlsx//;
 $file_merged .= $fii."_";
 foreach my $sheet (@{$excel -> {Worksheet}}) {
        $sheet -> {MaxRow} ||= $sheet -> {MinRow};
         foreach my $row (1 .. $sheet -> {MaxRow}) {
                my $cell_name = $sheet -> {Cells}[$row][0];
                my $cell_number = $sheet -> {Cells}[$row][1];
               if($cell_number){
#               	print $cell_name->{Val}." ".$cell_number->{Val}."\n";
               	$data{$sheet->{Name}}->{$fii}->{$cell_name->{Val}} = $cell_number->{Val};
                        }
                }
 }
}
#print Dumper(%data);
$file_merged .= ".xlsx";
my $workbook = Excel::Writer::XLSX -> new($file_merged);
my $bold      = $workbook->add_format( bold => 1 );

#my $scl = eval $data{'family_level'}->{'GZS1_2'}->{'Ruminococcaceae'};
#print $data{'family_level'}->{'GZS1_2_'}->{'Ruminococcaceae'}."\n";

my @levels = keys %data;
#print join(" ",@levels)."\n";

#    print Dump(%data);
foreach my $key(@levels){ ## HERE, if use $key(%data) will produce 2*(keys %data) $key, half real, half hash refs
#	print $key."\n";
my $worksheet = $workbook->add_worksheet($key);
    my @headings = keys %{$data{$key}};
#    print join(" ",@headings)."\n";
    my %all_names = ();
    foreach my $head(@headings){
    	foreach my $nam(sort keys %{$data{$key}->{$head}}){
    	$all_names{$nam} = 1;
    }
    }
#    print Dump{%all_names};
    my $count = 0;
    my @datadata = ();
    my @col1_name = ();
    my $nnn = $#headings+1;
    	foreach my $h(0..$nnn){
    	 if($h<$nnn){
        foreach my $nams(sort keys %all_names){
    		 unless(defined($data{$key}->{$headings[$h]}->{$nams})){ ## HERE, if use 'if(undef($data{$key}...))' will cause all number values lost
    			                                              ## cause undef is a function that undefines EXPR, undef can not be used for judgement
    			$data{$key}->{$headings[$h]}->{$nams} = 0;
    		 }
#    		 print $key." ".$head." ".$nams." ".$data{$key}->{$head}->{$nams}."\n";
    		 push @{$datadata[$count]},$data{$key}->{$headings[$h]}->{$nams};
    	}
    	$count++;
    }
    else{
    	foreach my $nams(sort keys %all_names){
    		push @col1_name,$nams;
    }
    }
  }
#    print Dump(%data);
#    print $data{'family_level'}->{'GZS1_2_'}->{'Ruminococcaceae'}."\n";
#    print Dump(@datadata);
    # Add the worksheet data that the charts will refer to.
    my $headings = [$key,@headings];
    my $datas = [
    [@col1_name],
        @datadata
    ];

    $worksheet->write( 'A1', $headings, $bold );
    $worksheet->write( 'A2', $datas );
    $worksheet->set_column(0,0,25);
}