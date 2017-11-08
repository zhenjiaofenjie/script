#! /usr/bin/perl

# This program can be used for sequence extraction.
# Based on extraction standard below:
# 1.Both up primers must be found. Up primers must be full-length.One single seq has no more than 1 up primer.
# 2.Barcodes are <barcode_length> bp length outside of the primers, they must be full-length. Up barcodes are required, further more you can choose to consider both barcodes in one sequence(be aware in this case downprimer must be full length too).
# 3.The length of VR region must be more than <min_length> bp and less than <max_length> bp.
# 4.VR region has no more than <N_number> N bases.
# 5.More than or equal to <good_percent> bases in one sequence must have better quality than <Qvalue>, or the sequence will be discarded.
# Default arguments are used to extract 454 16S rDNA V3 region sequences.
#
# Written by WangJing, wjingsjtu@gmail.com
#
# Changes
# v1.0
# 1.Change the steps, barcode check goes right behind primer check, so that it's exactly the same as extraction standard.
# 2.Change quality check from more than <good_percent> to more than and equal to <good_percent>.
# Modified in 2013.03.25
# v2.0
# 1.Add submit opt, if set, the original sequences belong to the extracted ones will be saved for submitting later
# 2.Change the default qvalue & good_percent to 0, means no quality check by default.
# Modified in 2013.04.07
# v3.0
# 1.Fix the %QualHash, so that can get the correct length of each qual sequence.
# Modified in 2013.06.04
# v1v3_v1.0
# 1.Change to fit v1v3 extraction.
# Modified in 2013.06.24


use strict;
use warnings;
use vars qw($USAGE);
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Bio::AlignIO;


$USAGE = "454_V3_wj.pl [arguments] --pf primerfile --sq sequencefile --ql qualityfile --tag tagfile
REQUIRED:
\t--pf/p\tprimerfile,include the forward and reverse primers, FASTA format
\t--sq/s\tsequencefile,FASTA format
\t--ql/q\tqualityfile,QUAL format
\t--tag/t\ttagfile,list sample id and barcode sequence, tab delimited
ARGUMENTS:
\t--h/help\tprint this text
\t--qn\tqiimenumber to start, each extracted sequence will be named as '<sample id>_<qiime number>', and qn increased by 1(DEFAULT:10000000)
\t--bl/barcode_length\tthe length of the barcode sequences(DEFAULT:6)
\t--min/min_length\tthe min length of the extracted region(DEFAULT:300)
\t--max/max_length\tthe max length of the extracted region(DEFAULT:700)
\t--N/N_number\tthe max number of the N bases in one sequence(DEFAULT:2)
\t--Qvalue\tthe Q-value threshold(DEFAULT:0)
\t--good/good_percent\tthe rate of the bases that should be better than Q-value threshold(DEFAULT:0)
\t--both\tuse to consider both up and down stream barcodes, sequences are extracted only when the two barcodes are exactly the same
\t--rp/revprimer\tthe name of reverse primer, if defined, the direction of seqs will be set to forward primer->seq region->reverse primer
\t--nosubmit\tunless it's set,or the sequences for submit, which without trimming of barcodes nor primers, will be saved for submitting to NCBI etc. later.\n";

my ($primerfile,$sequencefile,$tagfile,$qiimenum,$qualityfile,$help,$barcode_length,$min_length,$max_length,$N_number,$Qvalue,$good_percent,$both,$rev_primer,$submit) =
   (undef,      undef,        undef,   10000000,     undef,      "",   6,               300,       700,        2,        0,     0,           "",  undef,      1,);

&GetOptions('pf|p=s'             => \$primerfile,
	        'sq|s=s'             => \$sequencefile,
	        'tag|t=s'            => \$tagfile,
	        'qn=s'               => \$qiimenum,
	        'ql|q=s'             => \$qualityfile,
			'help|h'             => \$help,
            'barcode_length|bl=i'=> \$barcode_length,
            'min_length|min=i'   => \$min_length,
			'max_length|max=i'   => \$max_length,
			'N_number|N=i'       => \$N_number,
			'Qvalue=i'           => \$Qvalue,
			'good_percent|good=f'  => \$good_percent,
			'both'               => \$both,
			'revprimer|rp=s'     => \$rev_primer,
			'nosubmit!'			 => \$submit,
            );
if($help){
	die($USAGE);
}
if (!$primerfile || !$sequencefile || !$tagfile || !$qualityfile) {
    if(!$primerfile){
		print "Please specify a primer file!\n";
	}
	if(!$sequencefile){
		print "Please specify a sequence file!\n";
	}
	if(!$tagfile){
		print "Please specify a tag file!\n";
	}
	if(!$qualityfile){
		print "Please specify a qualityfile file!\n";
	}
	die("For more help please use --help/-h\n");
}

my $tm = localtime(time());
print "$tm\n454_extract_v1v3.pl\nThis is a beta version. If this program print some warnings, please contact me to see whether it's ok:)\n";
$tm=~ s/\s/_/g;
my $seq_out = $sequencefile;
$seq_out =~ s/fna|fasta|fa//;
$seq_out =$tm."_".$seq_out."fa";
my $submit_out=$seq_out;
$submit_out=~ s/\.fa/_submit.fa/;

my $log_out=$seq_out;
$log_out=~ s/\.fa/\.log/;

open(LOG,">".$log_out);   #write down the parameters used
print LOG "pf\t$primerfile\nsq\t$sequencefile\ntag\t$tagfile\nqn\tqiimenum\nql\t$qualityfile\nbarcode_length\t$barcode_length\nmin_length\t$min_length\nmax_length\t$max_length\nN_number\t$N_number\nQvalue\t$Qvalue\ngood_percent\t$good_percent\nboth\t$both\nrevprimer\t$rev_primer\nsubmit\t$submit";
close(LOG);

## Start of the program ##

## READ FASTA FILE
## Returns %SeqHash
my %SeqHash = ();
my ($key,$value) = ("","");
my $seq_count;
open FASTA, "$sequencefile" or die "Cannot open $sequencefile.\n";
while (<FASTA>)
{
	chomp;
	if (/>/)
	{
		$seq_count++;
		if ($key ne "" && $value ne "")
		{
			$SeqHash{$key} = $value;
			$key = "";
			$value = "";
		}
		my @qualary = split (" ", $_);
		$key = $qualary[0];
		$key =~ s/>//;
	}
	else
	{
		$value .=$_;
	}
}

if ($key ne "" && $value ne "")
{
	$SeqHash{$key} = $value;
	$key = "";
	$value = "";
}
close FASTA;
print "fasta Reading done.\n";

## READ QUALITY FILE
my %QualHash=();
($key,$value) = ("","");
open (QUAL,$qualityfile) or die "Cannot open $qualityfile.\n";
while(<QUAL>){
	chomp;
	if(/>/){
		if ($key ne "" && $value ne ""){
			$QualHash{$key}=$value;
			($key,$value)=("","");
		}
		my @line=split(" ",$_);
		$key=$line[0];
		$key=~ s/>//;
	}else{
		if($value eq ""){
			$value.=$_;
		}else{
			$value.=" ".$_; #Because there's a space between each position in qual file
		}
	}
}
if($key ne "" && $value ne ""){
	$QualHash{$key}=$value;
	($key,$value)=("","");
}
close(QUAL);
print "Quality file reading done.\n";

## READ TAG MAPPING FILE
my %tag_hash;
my $tag_count;
open TAGFILE,$tagfile or die"Can't Open tag file.\n";
foreach(<TAGFILE>){
	$tag_count++;
	chomp;
	my @tag_map = split;
	${$tag_hash{$tag_map[1]}}[0] =$tag_map[0];
	${$tag_hash{$tag_map[1]}}[1] =1;
}
close TAGFILE;

my $table_out = $tm."_".$sequencefile.".table.txt";
open TABLE,">$table_out";
print TABLE "RawSeqID\tSeqLength\tUpPrimer\tDnPrimer\tTooManyPrimers\tVRlength\tTag\tSample\tFinallyKeeped\tQualityNotFine\tBetterThan$Qvalue\tDiffUpDownTag\tUpTag\tDownTag\tNbasecount\tVR_avg_qual\n";
#    $res_print[0]       [1]       [2]      [3]        [4]               [5]   [6]     [7]    [8]            [9]             [10]               [11]          [12]   [13]	[14]       [15]
#Print the condition of each sequence into a table. The variable $res_print[] can be assigned to any important information.   
#The table contains 10 parts;

## Blastn the Primers
## formatdb function can also be added by system method
my $blastout = $tm.".out";
system("formatdb -i $primerfile -p F -o T/F");
my @para = (-database => $primerfile,-program => 'blastn', -e => 1e-1, -W => 4, -a => 4, -o =>$blastout);
my $factory = Bio::Tools::Run::StandAloneBlast->new(@para);
$factory->blastall($sequencefile);
print "Blast OK!\n";

## The main loop, process the blast results and produce extracted fasta file
open FOR, $blastout;
open(SEQOUT,">".$seq_out);
if($submit){
	open(SUBOUT,">".$submit_out);
}
my $report = Bio::SearchIO->new(-fh=>\*FOR, -format=>'blast');
my ($extracted_count,$upprimer_count,$dnprimer_count,$VRlength_ok_count,$up_tag_count,$down_tag_count,%sample_dis,$quality_ok_count,$tag_exist_count,$avg30_count);
while(my $result = $report->next_result){  # $result is a Bio::Search::Result::ResultI object
	my @res_print;
	my $revcomp=0;       #set to 1 if seq need to be revcomped
	$res_print[0]=$result->query_name();
	my $seqlength = $result->query_length();
	$res_print[1]=$seqlength;
	my $findup=0;
	my $finddown=0;
	my ($uptag,$uptag_pos,$downtag,$downtag_pos);
	my $VR_start_pos=0;
	my $VR_end_pos=$seqlength-1;  #position starts from 0 to the end
        while (my $hit = $result->next_hit()){ # # $hit is a Bio::Search::Hit::HitI object
		while( my $hsp = $hit->next_hsp ) {  # $hsp is a Bio::Search::HSP::HSPI object
			if($hsp->strand('hit')==1 && $hsp->length('hit')==$hit->length()){  #upprimer is found.  primers must be full length
				$findup++;
				$res_print[2]=$hsp->query_string;
				$VR_start_pos=$hsp->query->end;
				$uptag_pos=$hsp->query->start-$barcode_length-1;
				if($rev_primer && $rev_primer eq $hit->name()){   #if rev_primer is set, and upprimer is rev_primer, mark to revcomp this sequence later
					$revcomp=1;
				}
			}
			if($hsp->strand('hit')==-1 && !$both){              #downprimer is found.  primers may not be full length.
				$finddown++;
				$res_print[3]=$hsp->query_string;
				$VR_end_pos=$hsp->query->start-2;
			}
			if($hsp->strand('hit')==-1 && $both && $hsp->length('hit')==$hit->length()){              #downprimer is found.  primers must be full length.
				$finddown++;
				$res_print[3]=$hsp->query_string;
				$VR_end_pos=$hsp->query->start-2;
				$downtag_pos=$hsp->query->end;
			}
		}
	}
	if($findup>0){
		$upprimer_count++;
	}
	if($finddown>0){
		$dnprimer_count++;
	}
	unless($findup==1){     #upprimer must be found, and only one upprimer are found
		$res_print[4]=1;
		&table_print(@res_print);
		next;
	}
	$uptag=undef;
	if($uptag_pos>=0){       #tag must be full-length
		$uptag=substr($SeqHash{$result->query_name},$uptag_pos,$barcode_length);
		$res_print[12]=$uptag;
	}
	$downtag=undef;
	if($both){
		if($downtag_pos<=length($SeqHash{$result->query_name})-$barcode_length){        #tag must be full-length
			$downtag=substr($SeqHash{$result->query_name},$downtag_pos,$barcode_length);
			$downtag=reverse($downtag);   #the revcomp downtag is the real barcode
			$downtag=~ tr/AGTC/TCAG/;
			$res_print[13]=$downtag;
		}
		if(!$uptag || !$downtag || $uptag ne $downtag){            #when -both is set, up and down barcodes must be exactly the same
			$res_print[11]=1;
			&table_print(@res_print);			
			next;
		}
	}
	if(defined($uptag) && ${$tag_hash{$uptag}}[1]){
		$tag_exist_count++;
		$res_print[6]=$uptag;
		$res_print[7]=${$tag_hash{$uptag}}[0];
	}else{
		&table_print(@res_print);
		next;
	}
	my $VR_length=$VR_end_pos-$VR_start_pos+1;
	$res_print[5]=$VR_length;
	if($VR_length<$min_length || $VR_length>$max_length){      #length of VR region must be more than <min_length> and less than <max_length>
		&table_print(@res_print);
		next;
	}
	$VRlength_ok_count++;
	my $VR_seq=substr($SeqHash{$result->query_name},$VR_start_pos,$VR_length);
	if($revcomp==1){
		$VR_seq=reverse($VR_seq);
		$VR_seq=~ tr/ATGCatgcN/TACGtacgN/;
	}
	my $count_N=($VR_seq=~ tr/N/N/);
	$res_print[14]=$count_N;
	if($count_N>$N_number){       #number of N base in VR region must be no more than <N_number>
		&table_print(@res_print);
		next;
	}
	my @qual=split(/\s+/,$QualHash{$result->query_name});  #some how there're more than one spaces between positions
	my $bad_base=0;
	my $VR_qual_sum=0;
	for(my $i=0;$i<$seqlength;$i++){   #use $seqlength instead of @qual in case of the seqs in fasta having been trimed before
		if($qual[$i]<$Qvalue){
			$bad_base++;
		}
		if($i>=$VR_start_pos && $i<($VR_start_pos+$VR_length)){
			$VR_qual_sum+=$qual[$i];
		}
	}
	if((1-$bad_base/$seqlength)<$good_percent){          #More than or equal to <good_percent> bases in one sequence must have better quality than <Qvalue>
		$res_print[9]=1;
		$res_print[10]=1-$bad_base/$seqlength;
		&table_print(@res_print);
		next;
	}else{
		$quality_ok_count++;
	}
	$res_print[8]=1;
	$res_print[15]=$VR_qual_sum/$VR_length;
	if($VR_qual_sum/$VR_length>=30){
		$avg30_count++;
	}
	$sample_dis{${$tag_hash{$uptag}}[0]}++;
	print SEQOUT ">".${$tag_hash{$uptag}}[0]."_".$qiimenum." ".$result->query_name." orig_bc=".$uptag." new_bc=".$uptag." bc_differs=0 quality=\n".$VR_seq."\n";
	if($submit){
		print SUBOUT ">".${$tag_hash{$uptag}}[0]."_".$qiimenum." ".$result->query_name." orig_bc=".$uptag." new_bc=".$uptag." bc_differs=0 quality=\n".$SeqHash{$result->query_name}."\n";
	}
	$qiimenum++;
	$extracted_count++;
	&table_print(@res_print);
}
close FOR;
close SEQOUT;
close SUBOUT;

## Produce the report file
my $per_tag_extract=$extracted_count/$tag_count;
my $extract_per=$extracted_count/$seq_count*100;
my $upp_per=$upprimer_count/$seq_count*100;
my $dn_per=$dnprimer_count/$seq_count*100;
my $tag_OK_per=$tag_exist_count/$seq_count*100;
my $reppp = $seq_out;
my $VRlength_ok_per=$VRlength_ok_count/$seq_count*100;
my $quality_ok_per=$quality_ok_count/$seq_count*100;
$reppp =~ s/\.fa/\_report\.txt/;
open(REPORT,">".$reppp);

print REPORT "总序列数\t".$seq_count."\n样本数\t".$tag_count."\n提取出的序列数为\t".$extracted_count."\n平均提取量为\t".$per_tag_extract."\n提取出序列的比例(%)\t".$extract_per."\n找到上游引物的序列数\t".$upprimer_count."\n占序列总数比例(%)\t".$upp_per."\n找到下游引物的序列数\t".$dnprimer_count."\n占序列总数比例(%)\t".$dn_per."\ntag有效的序列数\t".$tag_exist_count."\n占序列总数比例(%)\t".$tag_OK_per."\n找到引物且可变区长度合适的序列数\t".$VRlength_ok_count."\n占序列总数比例(%)\t".$VRlength_ok_per."\n序列质量合适的序列数\t".$quality_ok_count."\n占序列总数比例(%)\t".$quality_ok_per."\n可变区平均质量大于30的序列数\t".$avg30_count."\n占序列总数比例(%)\t".($avg30_count/$seq_count*100)."\n";
close REPORT;

# Print sample distribution file
my $dis_out = $seq_out;
$dis_out =~ s/\.fa/\_sam_dis\.txt/;
open(DIS,">".$dis_out);
while((my $key,my $value)=each %sample_dis){
	print DIS "$key\t$value\n";
}
close DIS;
print STDOUT "All are done!\nThis is a beta version. If this program print some warnings, please contact me to see whether it's ok:)\n";
sub table_print{
	my @table_line = @_;
	for my $i(0..$#table_line){
		(defined($table_line[$i])) ? (print TABLE $table_line[$i]."\t") : (print TABLE "\t");
	}
	print TABLE "\n";
}
