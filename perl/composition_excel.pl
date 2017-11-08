#!/usr/bin/perl -w
## Pie graph showing microbiota-composition
## Written by Wang Guoyang
## paul.paulo.king@gmail.com
## 2011-11-03

## Problems to be solved
## 1.subclass and suborder, which level they belong to, how to show them on the graph, or just delete them?
## 2.colors of pie graph may not show right, must check before use
## 3.still not solve the rdp-online connection problem , also the stand-alone java program cannot work

## future features:
## solve above problems, and work with the extraction program
## thus form the one-step sample-analysis 'bigger program'

## Modified 2011-11-04
## subclass and suborder, etc, is removed from the original file, final pie will not show them


## BUGS #####

## Excel module cannot work well with Tk, maybe it's the mw competation
## also, Excel-relavent codes should be in the main part, not in the sub part

## cause some level have unclassified bacteria, so phylum-count summary may larger than class-count summary
## the unclassified bacteria is in phylum-count, not in class-count

## If bacteria name contains " ", may cause error, the reason is that when building a tree, it looks " " as level marks
## One Example is the 'Incertae Sedis XIII' and the 'Insertae Sedis XIV'

## Still, Notice TM7 phylum has only ONE genus, when the results show TM7, must pay attention to it

## ADVICE *****

## Double check is STRONGLY RECOMMANDED, as numbers-summary, some uncounted bacteria(with no numbers)
## and the excel should be re-formated to look better
## #########


## 2012-10-23
## When the Phylum and Class have the same name, eg. Actinobacteria
## The result xlsx can be strange


use strict;
use Excel::Writer::XLSX;
use Tree::Numbered::Tools;


my %hierarchy;
my $file = shift or die("USAGE: composition_excel.pl classification_file.txt\n");
my $nfile = "n_".$file;
my $fii = $file;
$fii =~ s/\.txt//;
my $leee_bef = 0;
my $leee;
my $tree;
my @phylum;
my @phylum_num;
my @class;
my @class_num;
my @order;
my @order_num;
my @family;
my @family_num;
my @genus;
my @genus_num;
my $workbook  = Excel::Writer::XLSX->new( $fii.'_composition.xlsx' );


readfile() or die "shoule have the classification file as parameters.\n";
reform_file();
form_tree();
select STDOUT;
pie_phylum();
pie_class();
pie_order();
pie_family();
pie_genus();


## Phylum chart
my $worksheet = $workbook->add_worksheet('phylum_level');
    my $bold      = $workbook->add_format( bold => 1 );

    # Add the worksheet data that the charts will refer to.
    my $headings = [ 'Phylum', 'Number' ];
    my $data = [
        [ @phylum ],
        [ @phylum_num     ],
    ];

    $worksheet->write( 'A1', $headings, $bold );
    $worksheet->write( 'A2', $data );
    $worksheet->set_column(0,0,25);

    # Create a new chart object. In this case an embedded chart.
    my $chart = $workbook->add_chart( type => 'pie', embedded => 1 );
    my $phy_n = @phylum;
#   $phy_n++;
    # Configure the series. Note the use of the array ref to define ranges:
    # [ $sheetname, $row_start, $row_end, $col_start, $col_end ].
    $chart->add_series(
        name       => 'Phylum '.$fii,
        categories => ['phylum_level',1,$phy_n,0,0],
        values     => ['phylum_level',1,$phy_n,1,1],
    );

    # Add a title.
    $chart->set_title( name => 'Phylum composition of '.$fii );

    # Set an Excel chart style. Colors with white outline and shadow.
    $chart->set_style( 10 );

    # Insert the chart into the worksheet (with an offset).
    $worksheet->insert_chart( 'C2', $chart, 8,6,1.8,1.5 );


## Class chart
$worksheet = $workbook->add_worksheet('class_level');

    # Add the worksheet data that the charts will refer to.
    $headings = [ 'Class', 'Number' ];
    $data = [
        [ @class ],
        [ @class_num     ],
    ];

    $worksheet->write( 'A1', $headings, $bold );
    $worksheet->write( 'A2', $data );
    $worksheet->set_column(0,0,25);

    # Create a new chart object. In this case an embedded chart.
    $chart = $workbook->add_chart( type => 'pie', embedded => 1 );
    my $cla_n = @class;
#   $phy_n++;
    # Configure the series. Note the use of the array ref to define ranges:
    # [ $sheetname, $row_start, $row_end, $col_start, $col_end ].
    $chart->add_series(
        name       => 'Class '.$fii,
        categories => ['class_level',1,$cla_n,0,0],
        values     => ['class_level',1,$cla_n,1,1],
    );

    # Add a title.
    $chart->set_title( name => 'Class composition of '.$fii );

    # Set an Excel chart style. Colors with white outline and shadow.
    $chart->set_style( 10 );

    # Insert the chart into the worksheet (with an offset).
    $worksheet->insert_chart( 'C2', $chart, 8,6,1.8,1.5 );


## Order chart
$worksheet = $workbook->add_worksheet('order_level');

    # Add the worksheet data that the charts will refer to.
    $headings = [ 'Order', 'Number' ];
    $data = [
        [ @order ],
        [ @order_num     ],
    ];

    $worksheet->write( 'A1', $headings, $bold );
    $worksheet->write( 'A2', $data );
    $worksheet->set_column(0,0,25);

    # Create a new chart object. In this case an embedded chart.
    $chart = $workbook->add_chart( type => 'pie', embedded => 1 );
    my $ord_n = @order;
#   $phy_n++;
    # Configure the series. Note the use of the array ref to define ranges:
    # [ $sheetname, $row_start, $row_end, $col_start, $col_end ].
    $chart->add_series(
        name       => 'Order '.$fii,
        categories => ['order_level',1,$ord_n,0,0],
        values     => ['order_level',1,$ord_n,1,1],
    );

    # Add a title.
    $chart->set_title( name => 'Order composition of '.$fii );

    # Set an Excel chart style. Colors with white outline and shadow.
    $chart->set_style( 10 );

    # Insert the chart into the worksheet (with an offset).
    $worksheet->insert_chart( 'C2', $chart, 8,6,1.8,1.5 );
    
## Family chart
$worksheet = $workbook->add_worksheet('family_level');

    # Add the worksheet data that the charts will refer to.
    $headings = [ 'Family', 'Number' ];
    $data = [
        [ @family ],
        [ @family_num     ],
    ];

    $worksheet->write( 'A1', $headings, $bold );
    $worksheet->write( 'A2', $data );
    $worksheet->set_column(0,0,25);

    # Create a new chart object. In this case an embedded chart.
    $chart = $workbook->add_chart( type => 'pie', embedded => 1 );
    my $fam_n = @family;
#   $phy_n++;
    # Configure the series. Note the use of the array ref to define ranges:
    # [ $sheetname, $row_start, $row_end, $col_start, $col_end ].
    $chart->add_series(
        name       => 'Family '.$fii,
        categories => ['family_level',1,$fam_n,0,0],
        values     => ['family_level',1,$fam_n,1,1],
    );

    # Add a title.
    $chart->set_title( name => 'Family composition of '.$fii );

    # Set an Excel chart style. Colors with white outline and shadow.
    $chart->set_style( 10 );

    # Insert the chart into the worksheet (with an offset).
    $worksheet->insert_chart( 'C2', $chart, 8,6,1.8,1.5 );


## Genus Chart
$worksheet = $workbook->add_worksheet('genus_level');

    # Add the worksheet data that the charts will refer to.
    $headings = [ 'Genus', 'Number' ];
    $data = [
        [ @genus ],
        [ @genus_num     ],
    ];

    $worksheet->write( 'A1', $headings, $bold );
    $worksheet->write( 'A2', $data );
    $worksheet->set_column(0,0,25);

    # Create a new chart object. In this case an embedded chart.
    $chart = $workbook->add_chart( type => 'pie', embedded => 1 );
    my $gen_n = @genus;
#   $phy_n++;
    # Configure the series. Note the use of the array ref to define ranges:
    # [ $sheetname, $row_start, $row_end, $col_start, $col_end ].
    $chart->add_series(
        name       => 'Genus '.$fii,
        categories => ['genus_level',1,$gen_n,0,0],
        values     => ['genus_level',1,$gen_n,1,1],
    );

    # Add a title.
    $chart->set_title( name => 'Genus composition of '.$fii );

    # Set an Excel chart style. Colors with white outline and shadow.
    $chart->set_style( 10 );

    # Insert the chart into the worksheet (with an offset).
    $worksheet->insert_chart( 'C2', $chart, 8,6,1.8,1.5 );


sub readfile{
	open HIE,$file;
}

sub reform_file{
	open OUT,">$nfile";
	select OUT;
	print "    level level_name number\n\n";
	foreach(<HIE>){
		chomp;
		my $line = $_;
		my $sub_count = 0;
		if($line =~ /([\S ]+?)(\w+[ \S]+)/){
			my $some = $1;
			my $here = $2;
			$leee = ($some=~tr/ / /);
			if($here =~ /(subclass)|(suborder)|(subfamily)/){
				$sub_count+=2;
				next;
			};
				my $sub_level = $leee_bef-$leee;
				if($sub_level>=0){ # means hierarchy goes to a ligher level
				my $sub_count -= $sub_level;
				if($sub_count<0){
				$sub_count = 0;}
			}
			$leee -= $sub_count;
			for (1..$leee){
				print " ";
			}
			print $here."\n";
			$leee_bef = $leee;
		}
}
close OUT;
}

sub form_tree{
	 $tree = Tree::Numbered::Tools->readFile(
              filename         => $nfile,
              usecolumnnames   => 1,
  );
}

sub pie_phylum{
	my %phylum;
#	my $tree = $tree->getSubTree(2);
	foreach($tree->listChildNumbers){
		my $no = $_;
		my @arr = $tree->follow($no,"level");
		my @names = $tree->follow($no,"level_name");
		my @numbers = $tree->follow($no,"number");
		my $arr = @arr;
		my $names = @names;
		my $numbers = @numbers;
		if($arr==2){
			if($arr[1] eq 'phylum'){
				my $numm = $numbers[1];
				$numm =~ /(\d+)/;
				$numm = $1;
			$phylum{$names[1]}=$numm;
		}
		else{
			my $numm = $names[1];
				$numm =~ /(\d+)/;
				$numm = $1;
			$phylum{$arr[1]}=$numm;
		 }
    }
  }
  foreach my $key(sort keys %phylum){
  	push @phylum,$key;
  	push @phylum_num,$phylum{$key};
  } 
 
}

sub pie_class{
#	my $root = $tree->getSubTree(1);
	my %class;
	foreach($tree->listChildNumbers){
		my $no = $_;
		my @arr = $tree->follow($no,"level");
		my @names = $tree->follow($no,"level_name");
		my @numbers = $tree->follow($no,"number");
		my $arr = @arr;
		if($arr==3){
			if($arr[2] eq 'class'){
				my $numm = $numbers[2];
				$numm =~ /(\d+)/;
				$numm = $1;
			$class{$names[2]}=$numm;
		}
		else{
				my $numm = $names[2];
				$numm =~ /(\d+)/;
				$numm = $1;
			$class{$arr[2]}=$numm;
		 }
    }
  }
  foreach my $key(sort keys %class){
  	push @class,$key;
  	push @class_num,$class{$key};
  }
    
}

sub pie_order{
#	my $root = $tree->getSubTree(1);
	my %order;
	foreach($tree->listChildNumbers){
		my $no = $_;
		my @arr = $tree->follow($no,"level");
		my @names = $tree->follow($no,"level_name");
		my @numbers = $tree->follow($no,"number");
		my $arr = @arr;
		if($arr==4){
			if($arr[3] eq 'order'){
				my $numm = $numbers[3];
				$numm =~ /(\d+)/;
				$numm = $1;
			$order{$names[3]}=$numm;
		}
		else{
				my $numm = $names[3];
				$numm =~ /(\d+)/;
				$numm = $1;
			$order{$arr[3]}=$numm;
		 }
    }
  }
  foreach my $key(sort keys %order){
  	push @order,$key;
  	push @order_num,$order{$key};
  }
}

sub pie_family{
#	my $root = $tree->getSubTree(1);
	my %family;
	foreach($tree->listChildNumbers){
		my $no = $_;
		my @arr = $tree->follow($no,"level");
		my @names = $tree->follow($no,"level_name");
		my @numbers = $tree->follow($no,"number");
		my $arr = @arr;
		if($arr==5){
			if($arr[4] eq 'family'){
				my $numm = $numbers[4];
				$numm =~ /(\d+)/;
				$numm = $1;
			$family{$names[4]}=$numm;
		}
		else{
				my $numm = $names[4];
				$numm =~ /(\d+)/;
				$numm = $1;
			$family{$arr[4]}=$numm;
		 }
    }
  }
  foreach my $key(sort keys %family){
  	if(defined($key)){
  	push @family,$key;
  push @family_num,$family{$key};}
  }
}

sub pie_genus{
#	my $root = $tree->getSubTree(1);
	my %genus;
	foreach($tree->listChildNumbers){
		my $no = $_;
		my @arr = $tree->follow($no,"level");
		my @names = $tree->follow($no,"level_name");
		my @numbers = $tree->follow($no,"number");
		my $arr = @arr;
		if($arr==6){
			if($arr[5] eq 'genus'){
				my $numm = $numbers[5];
				$numm =~ /(\d+)/;
				$numm = $1;
			$genus{$names[5]}=$numm;
		}
		else{
				my $numm = $names[5];
				$numm =~ /(\d+)/;
				$numm = $1;
			$genus{$arr[5]}=$numm;
		 }
    }
  }
  foreach my $key(sort keys %genus){
  	push @genus,$key;
  	push @genus_num,$genus{$key};
  }
    
}
