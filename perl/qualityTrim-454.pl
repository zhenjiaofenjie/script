#! /usr/bin/perl
use warnings;
$fas = shift;
$qual = shift;

$usage="USAGE: perl qualityTrim-454.pl fastafile qualfile\n";
my $sequence;
open (QUAL, "<$qual") or die($usage);
while (<QUAL>) {
    if($_ =~ /^>/) {
    	chomp;
	$name2 = $_;
	if($sequence ne "") {
	    @seq = split /\s+/, $sequence;
	    trim();
	}
	$sequence = "";
	$name = $name2;
    } else {
	$_ =~ s/\n$/ /;
	$sequence .= $_;
    }
}
@seq = split / /, $sequence;
trim();
close (QUAL);
print "quality trimmed.\n";


open OUT, ">trimed.fna";
select OUT;

$sequence = "";
open (FAS, "<$fas") or die;
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
    for(local $i=0; $i< @seq; $i++) { #check quality
	if($seq[$i] < 20) {
	    local $over = 0;# if a base has a quality less than 20, set over to 0
	    for(local $j= $i+1; $j < @seq; $j++) { # from the checked low quality base, check next bases
		if($j - $i >= 4) {
		    $over = 0;
		    last; # if the distance between j base and i base larger than 9, come out the loop, set over to 0
		    # means if there are 9 continuous low quality bases, thought it as end
		}
		if($seq[$j] >= 20) { # an example: i is 200, qual(i) is 19, then if j is 205, and qual(j) is 21, then j-i not bigger than 9, i is not end ,set i to 201
		    $over = 1;
		    last; # if j base has a quality higher than 20, come out the loop, set over to 1
		}
	    }
	    if($over == 0) {
		$end_cut = $i;
		last; # if over is 0, the i base is thought to be the end base, $end_cut set to i
	    }
	}
    }
    $remain{$name} = $end_cut; # returns a hash %remain, key is $name, value is $end_cut
}
