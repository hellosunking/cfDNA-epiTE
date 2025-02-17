#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <file.prefix> [target.size=4] [overall.motif.file]\n\n";
	exit 2;
}

my $min_freq = 1e-10;
my $norm = 1/log(256);
my $motifSize = $ARGV[1] || 4;

## guess sid
my $sid = $ARGV[0];
$sid =~ s/.*\///;

## get TE list
my @TE;
open IN, "$ARGV[0].bed.list" or die( "$!" );
while( <IN> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;
	push @TE, $l[0];
}
close IN;

## get diversity per TE
my @d;
foreach my $te ( @TE ) {
	if( -s "$ARGV[0].$te.motif" ) {
		open IN, "$ARGV[0].$te.motif" or die( "$!" );
	} elsif( -s "$te.motif" ) {
		open IN, "$te.motif" or die( "$!" );
	} else {
		push @d, 'NA';
		next;
	}

	my $v = 0;
	while( <IN> ) {
		chomp;
		my @l = split /\t/; ##CCCA	169229	5.51757888916132
		next unless length($l[0])==$motifSize && $l[0]=~/^[ACGT]+$/ && $l[-1] > $min_freq;
		$l[-1] /= 100;
		$v -= $l[-1] * log($l[-1]);
	}
	close IN;
	push @d, $v * $norm;
}

## overall? to be compatible with older versions
if( $#ARGV >= 2 ) {
	my $v = 0;
	open IN, "$ARGV[2]" or die( "$!" );
	while( <IN> ) {
		chomp;
		my @l = split /\t/;	##CCCA	395318	0.673892517963583	1582283	2.69729350801373
		next unless length($l[0])==$motifSize && $l[0]=~/^[ACGT]+$/ && $l[-1] > $min_freq;
		$l[-1] /= 100;
		$v -= $l[-1] * log($l[-1]);
	}
	close IN;
	#@m = map {$_ / $overall} @m;
	push @d, $v * $norm;
}
print join("\t", $sid, @d), "\n";


