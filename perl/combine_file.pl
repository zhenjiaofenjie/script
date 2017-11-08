#!/usr/bin/perl
#############################
#合并两个表格矩阵，行名相同的行合并，没有的部分用0补齐，最后去除和为0的行
#command line: perl combine_file.pl file1 file2 output_file
############################
# from Wang Jing
###########################
sub read_target_file{
	open(file,$ARGV[1]) or die("USAGE: combine_file.pl file1 file2 outfile\n");
	my @lines=<file>;
	close(file);
	chomp @lines;
	@length1=split(/\s+/,$lines[0]);
	$tittle=join("\t",@length1);
	shift(@length1);
	shift(@lines);
	$insert1=0;
	for(my $t=1;$t<@length1;$t++){
		$insert1=$insert1."\t0";
	}
	for(my $i=0;$i<@lines;$i++){
		my @col=split(/\s+/,$lines[$i]);
		my $j=$col[0];
		shift(@col);			
		$result{$j}=join("\t",@col);
	}
}
#############################

###########################
sub read_file{
	open(file,$ARGV[0]) or die("USAGE: combine_file.pl file1 file2 outfile\n");
	my @lines=<file>;
	close(file);
	chomp @lines;
	@length2=split(/\s+/,$lines[0]);
	shift(@length2);
	$tittle=$tittle."\t".join("\t",@length2);
	shift(@lines);
	$insert2=0;
	for(my $t=1;$t<@length2;$t++){
		$insert2=$insert2."\t0";
	}
	for(my $i=0;$i<@lines;$i++){
		my @col=split(/\s+/,$lines[$i]);
		my $j=$col[0];		
		if(!$result{$j}){
			$result{$j}=$result{$j}.$insert1;
		}
		shift(@col);
		$result{$j}=$result{$j}."\t".join("\t",@col);				
	}
}
#############################

##############################
sub write_file{
	system("touch $ARGV[2]") or die("USAGE: combine_file.pl file1 file2 outfile\n");	
	open(file,">".$ARGV[2]);
	print file $tittle."\n";
	foreach $uid(keys(%result)){
		my @a=split(/\s+/,$result{$uid});
		my $sum=0;
		for(my $count=0;$count<@a;$count++){
			$sum=$sum+$a[$count];
		}
		if($sum!=0){
			if(@length1+@length2>@a){
				print file $uid."\t".$result{$uid}."\t".$insert2."\n";			
			}else{
				print file $uid."\t".$result{$uid}."\n";
			}			
		}
	}
	close(file);	
}
#############################

read_target_file();
read_file();
write_file();
