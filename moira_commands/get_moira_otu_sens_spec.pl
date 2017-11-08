#!/usr/bin/perl
use strict;
use warnings;

###################
#align filtered seqs to mock reference, find out which seqs are good
###################
#system('~/mothur/mothur "#align.seqs(fasta=filtered.resize.fasta,reference=~/Illumina_Miseqrun/run04_20141113/run03mock_reference_21.align,processors=12)"');
open(SEQREPORT,"filtered.resize.align.report");
my %seqgood;
<SEQREPORT>;
while(<SEQREPORT>){
    chomp;
    my @line=split(/\t/,$_);
    my $seqname=shift @line;
    my $identity=pop @line;
    if($identity>=97){
        $seqgood{$seqname}=1;
    }
}
close(SEQREPORT);

open(OTUSENSSPEC,">otu_sens_spec.txt");
print OTUSENSSPEC "Method\tNo.ofOTUs\tTpmock\tTnmock\tFpmock\tFnmock\tTp\tTn\tFp\tFn\tMCC\n";

system('~/mothur/mothur "#align.seqs(fasta=usearch/mock_usearch_pcrstriped_otus_numbered.fasta,reference=~/Illumina_Miseqrun/run04_20141113/run03mock_reference_21.align,processors=12)"');
system('~/mothur/mothur "#align.seqs(fasta=qiime/rep_set/filtered.resize.nogap_rep_set.fasta,reference=~/Illumina_Miseqrun/run04_20141113/run03mock_reference_21.align,processors=12)"');
system('~/mothur/mothur "#align.seqs(fasta=mothur/filtered.resize.precluster.an.unique_list.0.03.rep.fasta,reference=~/Illumina_Miseqrun/run04_20141113/run03mock_reference_21.align,processors=12)"');
system("cp usearch/mock_usearch_pcrstriped_otus_numbered.align.report usearch/usearch.align.report");
system("cp qiime/rep_set/filtered.resize.nogap_rep_set.align.report qiime/qiime.align.report");
system("cp mothur/filtered.resize.precluster.an.unique_list.0.03.rep.align.report mothur/mothur.align.report");

###################
#read count table of unique seqs
###################
my %seqcount;
open(COUNT,"filtered.resize.count_table");
while(<COUNT>){
    chomp;
    my @line=split(/\t/,$_);
    my $seqname=shift @line;
    shift @line;
    @{$seqcount{$seqname}}=@line;
}
close(COUNT);
###################
#calculate how many OTUs and how many seqs can be remapped to OTUs
###################
my @remapnum;
my @otuhashes;
my @otunumbers;
my @methodnames=("usearch","qiime","mothur");

####################
#reading usearch readmap.uc
####################
open(UC,"usearch/mock_usearch_pcrstriped_readmap.uc");
my %otuassign;
while(<UC>){
    chomp;
    if(/^H/){
        my @line=split(/\t/,$_);
        push(@{$otuassign{$line[9]}},$line[8]);
    }
}
push(@otuhashes,{%otuassign});
close(UC);
###################
#make list file
###################
my $listname="usearch/usearch.list";
open(LIST,">".$listname);
my $otunum=keys(%otuassign);
push(@otunumbers,$otunum);
print LIST "0.03\t$otunum";
while(my($key,$value)=each(%otuassign)){
    print LIST "\t".join(",",@{$value});
}
close(LIST);

####################
#make usearch OTU table
####################
open(OTU,">usearch/otu_table.txt");
print OTU "Representative_Sequence\t".join("\t",@{$seqcount{'Representative_Sequence'}})."\n";
while((my $key,my $value)=each(%otuassign)){
    my @abund;
    foreach my $seqname(@{$value}){
        for(my $i=0;$i<@{$seqcount{$seqname}};$i++){
            $abund[$i]+=${$seqcount{$seqname}}[$i];
        }
    }
    print OTU $key."\t".join("\t",@abund)."\n";
}
close(OTU);
system('biom convert -i usearch/otu_table.txt -o usearch/otu_table.biom --to-json');

####################
#reading qiime uclust.otus.txt
####################
%otuassign=();
open(UC,"qiime/uclust_picked_otus/filtered.resize.nogap_otus.txt");
while(<UC>){
    chomp;
    my @line=split(/\t/,$_);
    my $otuname=shift @line;
    @{$otuassign{$otuname}}=@line;
}
push(@otuhashes,{%otuassign});
close(UC);
###################
#make list file
###################
$listname="qiime/qiime.list";
open(LIST,">".$listname);
$otunum=keys(%otuassign);
push(@otunumbers,$otunum);
print LIST "0.03\t$otunum";
while(my($key,$value)=each(%otuassign)){
    print LIST "\t".join(",",@{$value});
}
close(LIST);
####################
#make qiime OTU table
####################
open(OTU,">qiime/otu_table.txt");
print OTU "Representative_Sequence\t".join("\t",@{$seqcount{'Representative_Sequence'}})."\n";
while((my $key,my $value)=each(%otuassign)){
    my @abund;
    foreach my $seqname(@{$value}){
        for(my $i=0;$i<@{$seqcount{$seqname}};$i++){
            $abund[$i]+=${$seqcount{$seqname}}[$i];
        }
    }
    print OTU $key."\t".join("\t",@abund)."\n";
}
close(OTU);
system('biom convert -i qiime/otu_table.txt -o qiime/otu_table.biom --to-json');

####################
#reading mothur list file
####################
system('~/mothur/mothur "#align.seqs(fasta=mothur/filtered.resize.precluster.fasta,reference=~/Illumina_Miseqrun/run04_20141113/run03mock_reference_21.align,processors=12)"');
open(REP,"mothur/filtered.resize.precluster.an.unique_list.0.03.rep.fasta");
my %otu2rep;
while(<REP>){
    if(/^>/){
        my @line=split(/\t/,$_);
        my @name=split(/\|/,$line[1]);
        $line[0]=~ s/^>//;
        $otu2rep{$name[0]}=$line[0];
    }
}
close(REP);

open(LIST,"mothur/filtered.resize.precluster.an.unique_list.list");
%otuassign=();
my @otuname;
while(<LIST>){
    chomp;
    my @line=split(/\t/,$_);
    if(/numOtus/){
        shift @line;
        shift @line;
        @otuname=@line;
    }elsif($line[0] eq 0.03){
        shift @line;
        shift @line;
        foreach my $seqname(@line){
            @{$otuassign{$otu2rep{shift(@otuname)}}}=split(/,/,$seqname);
        }
    }
}
$otunum=keys(%otuassign);
push(@otunumbers,$otunum);
push(@otuhashes,{%otuassign});
close(LIST);

####################
#make mothur otu table
####################
system('~/mothur/mothur "#make.biom(shared=mothur/filtered.resize.precluster.an.unique_list.shared)"');

###################
#based on the identity to reference 
#TP:seq>=97% & OTU>=97%
#TN:seq<97% & OTU<97%
#FP:seq<97% & OTU>=97%
#FN:seq>=97% & OTU<97%
###################
my(@mocktps,@mocktns,@mockfps,@mockfns);
foreach my $method(@methodnames){
    my %otuassign=%{shift @otuhashes};
    my($tp,$tn,$fp,$fn)=(0,0,0,0);
    open(OTUREPORT,$method."/".$method.".align.report");
    <OTUREPORT>;
    while(<OTUREPORT>){
        chomp;
        my @line=split(/\t/,$_);
        my $otuname=shift @line;
        my $identity=pop @line;
        if(!$otuassign{$otuname}){
            next;
        }
        if($identity>=97){
            my @seqnames=@{$otuassign{$otuname}};
            foreach my $seqname(@seqnames){
                if($seqgood{$seqname}){
                    $tp++;
                }else{
                    $fp++;
                }
            }
        }else{
            my @seqnames=@{$otuassign{$otuname}};
            foreach my $seqname(@seqnames){
                if($seqgood{$seqname}){
                    $fn++;
                }else{
                    $tn++;
                }
            }
        }
    }
    close(OTUREPORT);
    push(@mocktps,$tp);
    push(@mocktns,$tn);
    push(@mockfps,$fp);
    push(@mockfns,$fn);
}

###################
#calculate MCC based on the pairwise distance of seqs
#TP:similarity between seqs>=97% & in the same OTU
#TN:similarity between seqs<97% & in different OTUs
#FP:similarity between seqs<97% & in the same OTU
#FN:similarity between seqs>=97% & in different OTUs
###################
my (@tps,@tns,@fps,@fns,@mccs);
open(PARALLEL,">para.txt");
print PARALLEL "usearch/usearch.list\tfiltered.resize.dist\nqiime/qiime.list\tfiltered.resize.dist\nmothur/filtered.resize.precluster.an.unique_list.list\tmothur/filtered.resize.precluster.dist";
close(PARALLEL);
system("parallel --gnu -a para.txt -C '\t' sens_spec.pl {1} {2} {1.}.sensspec");
open(SENSSPEC,"usearch/usearch.sensspec");
<SENSSPEC>;
$_=<SENSSPEC>;
chomp $_;
my @line=split(/\t/,$_);
push(@tps,$line[0]);
push(@tns,$line[1]);
push(@fps,$line[2]);
push(@fns,$line[3]);
push(@mccs,$line[4]);
close(SENSSPEC);
open(SENSSPEC,"qiime/qiime.sensspec");
<SENSSPEC>;
$_=<SENSSPEC>;
chomp $_;
@line=split(/\t/,$_);
push(@tps,$line[0]);
push(@tns,$line[1]);
push(@fps,$line[2]);
push(@fns,$line[3]);
push(@mccs,$line[4]);
close(SENSSPEC);
open(SENSSPEC,"mothur/filtered.resize.precluster.an.unique_list.sensspec");
<SENSSPEC>;
$_=<SENSSPEC>;
chomp $_;
@line=split(/\t/,$_);
push(@tps,$line[0]);
push(@tns,$line[1]);
push(@fps,$line[2]);
push(@fns,$line[3]);
push(@mccs,$line[4]);
close(SENSSPEC);
###################
#print the *_threshold_otu_num.txt, including num of OTUs, seqs being remmaped, MCCmock and MCCdist
###################
$otunum=shift @otunumbers;
my $remap=shift @remapnum;
my $tpmock=shift @mocktps;
my $tnmock=shift @mocktns;
my $fpmock=shift @mockfps;
my $fnmock=shift @mockfns;
my $tp=shift @tps;
my $tn=shift @tns;
my $fp=shift @fps;
my $fn=shift @fns;
my $mcc=shift @mccs;
print OTUSENSSPEC "USEARCH\t$otunum\t$tpmock\t$tnmock\t$fpmock\t$fnmock\t$tp\t$tn\t$fp\t$fn\t$mcc\n";
$otunum=shift @otunumbers;
$remap=shift @remapnum;
$tpmock=shift @mocktps;
$tnmock=shift @mocktns;
$fpmock=shift @mockfps;
$fnmock=shift @mockfns;
$tp=shift @tps;
$tn=shift @tns;
$fp=shift @fps;
$fn=shift @fns;
$mcc=shift @mccs;
print OTUSENSSPEC "QIIME\t$otunum\t$tpmock\t$tnmock\t$fpmock\t$fnmock\t$tp\t$tn\t$fp\t$fn\t$mcc\n";
$otunum=shift @otunumbers;
$remap=shift @remapnum;
$tpmock=shift @mocktps;
$tnmock=shift @mocktns;
$fpmock=shift @mockfps;
$fnmock=shift @mockfns;
$tp=shift @tps;
$tn=shift @tns;
$fp=shift @fps;
$fn=shift @fns;
$mcc=shift @mccs;
print OTUSENSSPEC "mothur\t$otunum\t$tpmock\t$tnmock\t$fpmock\t$fnmock\t$tp\t$tn\t$fp\t$fn\t$mcc\n";
chdir("..");
