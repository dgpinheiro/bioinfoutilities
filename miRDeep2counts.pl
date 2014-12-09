#!/usr/bin/perl

use strict;
use warnings;

my $infile = $ARGV[0];

die "Missing input file" unless ($infile);
die "Wrong input file ($infile)" unless (-e $infile);

open(IN, "<", $infile) or die "$!";

my $header_line = <IN>;
chomp($header_line);

my @header = split(/\t/, $header_line);

my @sample;
foreach my $h (@header) {
	if ($h !~ /^(#miRNA|read_count|precursor|total|.*\(norm\))$/) {
        $h=~s/^ *//g;
        $h=~s/ *$//g;
		push(@sample, $h);
	}
}

my %counts;

while(<IN>) {
	chomp;
	my %data;
	@data{ @header } = split(/\t/, $_);
	foreach my $s (@sample) {
		$counts{ $data{ $header[0] } }->{ $s } = 0 if ((! exists $counts{ $data{ $header[0] } })&&(! exists $counts{ $data{ $header[0] } }->{ $s }));
		$counts{ $data{ $header[0] } }->{ $s }+=$data{ $s };
	}
}

close(IN);

print join("\t", 'mirna', @sample),"\n";
foreach my $k (keys %counts) {
	print join("\t", $k, @{ $counts{ $k } }{ @sample }),"\n";
}
