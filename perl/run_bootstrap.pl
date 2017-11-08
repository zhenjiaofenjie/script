#!/usr/bin/perl
use strict;
use warnings;

for(my $i=0;$i<100;$i++){
	system("~/sparcc/SparCC.py boot/boot_$i.txt -i 100 -c boot/sim_cor_$i.txt >>sparcc.log &");
}
