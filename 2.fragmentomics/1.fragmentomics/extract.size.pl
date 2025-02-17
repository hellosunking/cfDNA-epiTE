#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <file.prefix> [target=150] [overall.size.file]\n\n";
	exit 2;
}

my $target = $ARGV[1] || 150;

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

## get size
my @s;
foreach my $te ( @TE ) {
	if( -s "$ARGV[0].$te.size" ) {
		open IN, "$ARGV[0].$te.size" or die( "$!" );
	} elsif( -s "$te.size" ) {
		open IN, "$te.size" or die( "$!" );
	} else {
		push @s, 'NA';
		next;
	}

	while( <IN> ) {
		chomp;
		my @l = split /\t/; ##150	70422	0.945116488808096	0.178988329981499	extra
		if( $l[0] eq $target ) {
			push @s, $l[3]*100;
			last;
		}
	}
	close IN;
}

## overall? to be compatible with older versions
if( $#ARGV >= 2 ) {
	my $overall = 0;
	open IN, "$ARGV[2]" or die( "$!" );
	while( <IN> ) {
		chomp;
		my @l = split /\t/;	##150	633037	0.937966520315893	0.180104708876802
		$overall = $l[-1]*100 if $l[0] eq $target;
	}
	close IN;
	#@s = map {$_ / $overall} @s;
	push @s, $overall;
}
print join("\t", $sid, @s), "\n";

