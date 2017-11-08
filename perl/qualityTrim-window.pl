#! /usr/bin/perl
use Getopt::Long;
use warnings;

my $fas = shift;
my $qual = shift;


my $outfile=$fas;
$outfile=~ s/(\.fasta|\.fna|\.fa)/_trimed\.fasta/;
my($window_size,$Qvalue)=
(  50,          20);

&GetOptions('windowsize|w=s'             => \$window_size,
	        'Qvalue|Q=s'             => \$Qvalue,
            );

$usage="USAGE: perl qualityTrim-window.pl fastafile qualfile
ARGUMENTS:
\t--windowsize/w\tEnable sliding window test of quality scores. If the 
                        average score of a continuous set of w nucleotides
                        falls below the threshold (see -Q for default), the
                        sequence is discarded. A good value would be 50. 
                        Default behavior for this function is to truncate the
                        sequence at the beginning of the poor quality window[default: 50]
\t--Qvalue/Q\tmin average qual score allowed in read [default: 25]
";

open (QUAL,$qual) or die($usage);
while (<QUAL>) {
    if($_ =~ /^>/) {
    	chomp;
	$name2 = $_;
	if($sequence && $sequence ne "") {
	    @seq = split /\s+/, $sequence;
	    trim();
	}
	$sequence = "";
	$name = $name2;
    } else {
	$_ =~ s/\n/ /;
	$sequence .= $_;
    }
}
@seq = split / /, $sequence;
trim();
close (QUAL);
print "quality trimmed.\n";


open OUT, ">".$outfile;
select OUT;
$sequence = "";
open (FAS,$fas) or die($usage);
while (<FAS>) {
    if($_ =~ /^>/) {
    	chomp;
	$name2 = $_;
	if($sequence ne "") {
	    $sequence = substr($sequence, 0, $remain{$name});
	    print "$name\n$sequence\n";
	}
	$sequence = "";
	$name = $name2;
    } else {
	chomp $_;
	$sequence .= $_;
    }
}

$sequence = substr($sequence, 0, $remain{$name}); # here remains these sequences that has a better 'end'
print "$name\n$sequence\n";


sub trim {
    local $end_cut = scalar(@seq); #$end_cut is seq length
	local $sumqual=0;
	local $averagequal=0;
	for(local $i=0; $i+$window_size<= @seq; $i++) { #check quality
		if($i==0){
			for(local $j= 0; $j < $window_size; $j++) { #calculate the average quality of first window start from first base.
				$sumqual+=$seq[$j];
			}	
		}else{
			$sumqual=$sumqual-$seq[$i-1]+$seq[$i+$window_size-1]; #window slide 1 step to right.
		}
		$averagequal=$sumqual/$window_size;
	    if($averagequal <$Qvalue) {
		$end_cut = $i;
		last; # if average quality of this window is less than threshold, $i base is thought to be the end base, $end_cut set to $i
	    }
    }
    $remain{$name} = $end_cut; # returns a hash %remain, key is $name, value is $end_cut
}
