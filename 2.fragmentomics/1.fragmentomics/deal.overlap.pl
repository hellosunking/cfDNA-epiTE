#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <cfDNA.ol.Repeat.wao>\n\n";
	exit 2;
}

my $minQ = 30;

my $lastLoci = -1;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr1	10008	10109	0	-	.	-1	-1	0
	next if $l[3] < $minQ || $l[5] eq '.';

	next if $l[1] == $lastLoci;
	## keep reads with ENDS located in the repeat region
	if( ($l[1]>=$l[6] && $l[1]<$l[7]) || ($l[2]>$l[6] && $l[2]<=$l[7]) ) {
		print "$l[0]\t$l[1]\t$l[2]\n";
		$lastLoci = $l[1];
	}
}
close IN;

