#!/usr/bin/perl

########################
# Daniel Guariz Pinheiro
# dgpinheiro@gmail.com


use strict;
use warnings;

my $samfile=$ARGV[0];

die "Missing SAM file" unless ($samfile);
die "Wrong SAM file" unless (-e $samfile);

open(SAM, "<", $samfile) or die $!;

my %data;
while(<SAM>) {
	chomp;
	if ($_=~/^@/) {
		print $_,"\n";
	} else {
		my @d = split(/\t/, $_);
		push(@{ $data{$d[0]} }, \@d);
	}
}

close(SAM);

foreach my $name (keys %data) {
	if (scalar( @{ $data{$name} }) == 1) {
		print join("\t", @{ $data{$name}->[0] } ),"\n";
	}
}
