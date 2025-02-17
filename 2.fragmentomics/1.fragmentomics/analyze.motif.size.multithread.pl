#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL
#

use strict;
use warnings;
use threads;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <genome.fa> <file.list> [thread=8]\n";
	print STDERR "\nCalculate 5'-NNNN 1-4mer end motifs.\n\n";
	exit 2;
}

my $threadNum = $ARGV[2] || 8;

my $autosomeOnly = 1;	## enforced!!!
my $min_mapQ  : shared = 30;
my $ACGT_only : shared = 1;	## only consider A/C/G/T, discard all other letters (N,Y,W,etc.)

print "Loading genome $ARGV[0] ...\n";
my $g : shared = load_genome( $ARGV[0] );

my (@sids, @fileLists);
open LST, "$ARGV[1]" or die( "$!" );
while( <LST> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	# could be '/path/to/bed' or 'SID /path/to/bed'
	push @sids, $l[0];
	push @fileLists, $l[-1];
}
close LST;

my $fileCnt = $#fileLists + 1;
for( my $i=0; $i<$fileCnt; $i+=$threadNum ) {
	my $j = $i + $threadNum;
	$j = $fileCnt if $j > $fileCnt;
	print STDERR "Batch: $i => $j\n";

	my @workers;
	for( my $k=$i; $k!=$j; ++$k ) {
		#print STDERR "$k: $fileLists[$k]\n";
		push @workers, threads->new(\&analyzer, $sids[$k], $fileLists[$k]);
	}

	foreach my $w ( @workers ) {
		$w->join;
	}
}

sub analyzer {
	my $sid   = shift;
	my $ifile = shift;

#	print "\rLoading $ifile ...\n";
	if( $ifile =~ /\.gz$/ ) {
		open IN, "zcat $ifile |" or die( "$!: $ifile" );
	} else {
		open IN, "$ifile" or die( "$!: $ifile" );
	}

	my $allcnt = 0;
	my %size;
	my (%m1, %m2, %m3, %m4);
	my ($c1, $c2, $c3, $c4) = ( 0, 0, 0, 0 );
	my ($s1, $s2, $s3, $s4);
	my %s2m;
	while( <IN> ) {
		chomp;
		my @l = split /\t/;

		next if $#l>=3 && $l[3]<$min_mapQ;
		next unless $l[0]=~/\d$/;	## autosome only
		next unless exists $g->{$l[0]};
		$s4 = substr( $g->{$l[0]}, $l[1], 4 );

		if( $ACGT_only ) {	## check N
			next unless $s4 =~ /^[ACGT]+$/;
		}
		
		my $here = $l[2] - $l[1];
		++ $size{$here};
		++ $allcnt;

		my $s3 = substr( $s4, 0, 3 );
		my $s2 = substr( $s4, 0, 2 );
		my $s1 = substr( $s4, 0, 1 );

		++ $c1;++ $m1{$s1};
		++ $c2;++ $m2{$s2};
		++ $c3;++ $m3{$s3};
		++ $c4;++ $m4{$s4};

		++ $s2m{"$here:$s1"};
	}
	close IN;

	my $ofile = $sid;
	$ofile =~ s/^.*\///;
	$ofile =~ s/\.gz$//;
	$ofile =~ s/\.bed$//;
	open OUT, ">$ofile.motif" or die( "$!" );
	print OUT "#Motif\tCount\tFrequency%\n";
	foreach my $i ( sort keys %m1 ) {
		print OUT join("\t", $i, $m1{$i}, $m1{$i}/$c1*100), "\n";
	}
	foreach my $i ( sort keys %m2 ) {
		print OUT join("\t", $i, $m2{$i}, $m2{$i}/$c2*100), "\n";
	}
	foreach my $i ( sort keys %m3 ) {
		print OUT join("\t", $i, $m3{$i}, $m3{$i}/$c3*100), "\n";
	}
	foreach my $i ( sort keys %m4 ) {
		print OUT join("\t", $i, $m4{$i}, $m4{$i}/$c4*100), "\n";
	}
	close OUT;

	open OUT, ">$ofile.size" or die( "$!" );
	print OUT "#Size\tCount\tFrequency%\tCumulative\tA%\tC%\tG%\tT%\n";
	my $cumu = 0;
	foreach my $i ( sort {$a<=>$b} keys %size ) {
		my $here = $size{$i} / $allcnt;
		$cumu += $here;
		my ($a, $c, $g, $t) = ( 0,0,0,0 );
		$a = $s2m{"$i:A"}/$size{$i}*100 if exists $s2m{"$i:A"};
		$c = $s2m{"$i:C"}/$size{$i}*100 if exists $s2m{"$i:C"};
		$g = $s2m{"$i:G"}/$size{$i}*100 if exists $s2m{"$i:G"};
		$t = $s2m{"$i:T"}/$size{$i}*100 if exists $s2m{"$i:T"};
		print OUT join("\t", $i, $size{$i}, $here*100, $cumu, $a, $c, $g, $t), "\n";
	}
	close OUT;
}

sub load_genome {
	my $fasta = shift;

	my %g;
	my $chr = 'NULL';
	open IN, "$fasta" or die("$!");
	while( <IN> ) {
		chomp;
		if( /^>(\S+)/ ) {
			$chr = $1;
			$g{$chr} = "";	## placeholder
		} else {
			$g{$chr} .= uc $_;
		}
	}
	close IN;
	return \%g;
}

