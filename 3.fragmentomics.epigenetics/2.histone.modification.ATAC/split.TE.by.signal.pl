#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date: 241220, consider the size of TEs for normalization

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.TE.ol.epigenome.cnt> <out.prefix>\n\n";
	exit 2;
}

my $sizeNormFactor = 1000;

my %TE2cnt;
my @normCnt;

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	## chr start end extra cnt
	my $v = $l[-1] / ($l[2]-$l[1]) * $sizeNormFactor;
	push @normCnt, $v;

	my $key = "$l[0]\t$l[1]\t$l[2]";
	$TE2cnt{$key} = $v;
}
my $total = $.;
close IN;

my @sortedCnt = sort {$a <=> $b} @normCnt;
my $cutoff = $sortedCnt[ $#normCnt >> 1 ];

open LOG, ">$ARGV[1].split.log" or die( "$!" );
print LOG "File\t$ARGV[0]\nTotal\t$total\nCutoff\t$cutoff\n";

open L, "| sort -k1,1 -k2,2n >$ARGV[1].low.bed" or die( "$!" );
open H, "| sort -k1,1 -k2,2n >$ARGV[1].high.bed" or die( "$!" );

my ($nLow, $nHigh) = (0, 0);
foreach my $k ( keys %TE2cnt ) {
	if( $TE2cnt{$k} <= $cutoff ) {
		++ $nLow;
		print L "$k\n";
	} else {
		++ $nHigh;
		print H "$k\n";
	}
}
close L;
close H;

print LOG "Low\t$nLow\nHigh\t$nHigh\n";
close LOG;
