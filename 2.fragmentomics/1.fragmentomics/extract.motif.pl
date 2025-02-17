#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <file.prefix> [target=CCCA] [overall.motif.file]\n\n";
	exit 2;
}

my $target = $ARGV[1] || 'CCCA';
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

## get motif per TE
my @m;
foreach my $te ( @TE ) {
	if( -s "$ARGV[0].$te.motif" ) {
		open IN, "$ARGV[0].$te.motif" or die( "$!" );
	} elsif( -s "$te.motif" ) {
		open IN, "$te.motif" or die( "$!" );
	} else {
		push @m, 'NA';
		next;
	}

	while( <IN> ) {
		chomp;
		my @l = split /\t/; ##CCCA	169229	5.51757888916132
		if( $l[0] eq $target ) {
			push @m, $l[-1];
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
		my @l = split /\t/;	##CCCA	395318	0.673892517963583	1582283	2.69729350801373
		$overall = $l[-1] if $l[0] eq $target;
	}
	close IN;
	#@m = map {$_ / $overall} @m;
	push @m, $overall;
}
print join("\t", $sid, @m), "\n";


