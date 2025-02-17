#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <TE.info> <genome.info> <file.prefix> <sampleid>\n";
	print STDERR "\nThe coverages are normalized to element length and sequencing depth.\n\n";
	exit 2;
}

## guess sid
my $sid = $ARGV[2];
$sid =~ s/.*\///;

## get TE list
my @TE;
open IN, "$ARGV[2].bed.list" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;
	push @TE, $l[0];
}
close IN;

## size of TE
my %TEsize;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##DNA	hAT-Charlie	268067	43492388
	$TEsize{"$l[0]"} = $l[-1];
}
close IN;

## overall
my $genomeSize = 0;
open IN, "$ARGV[1]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/; ##chr len
	$genomeSize += $l[-1] if $l[0]=~/^chr\d+$/;	## autosome only
}
close IN;
$TEsize{"Overall"} = $genomeSize;

## TE coverage
my @c;
foreach my $te ( @TE ) {
	open IN, "$te.size" or die( "$!" );
	<IN>;	## header
	my $read = 0;
	while( <IN> ) {
		chomp;
		my @l = split /\t/; ##150	70422	0.945116488808096	0.178988329981499	extra
		$read += $l[1];
	}
	close IN;
	$te =~ s/^$ARGV[3]\.//;
	push @c, $read/$TEsize{$te}*10000;
}

## normalize to overall, which is at the tail of @c
my $overall = pop @c;
@c = map { $_/ $overall } @c;

print join("\t", $sid, @c), "\n";

