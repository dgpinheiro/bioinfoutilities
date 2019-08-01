#!/usr/bin/env perl

use strict;
use warnings;

my $fastqfile = $ARGV[0];

die "Missing single-end .fastq file" unless ($fastqfile);
die "Wrong single-end .fastq file ($fastqfile)" unless (-e $fastqfile);

my $id;
my $desc;
my $seq;
my $qual;

open(IN, "<", $fastqfile) or die $!;
while(<IN>) {
    chomp;
    if ($. % 4 == 1) {
        ($id, $desc)=$_=~/^(\S+)(?:\s*(.*))?/;
    } elsif ($. % 4 == 2) {
        $seq=$_;
    } elsif ($. % 4 == 0) {
        $qual=$_;

        my $newqual=$qual;
        my $newseq=$seq;
        
        print   $id,(($desc) ? " $desc" : ''),"\n",
                $newseq,"\n",
                '+',"\n",
                $newqual,"\n";

        my $newqualrev = reverse($newqual);
        my $newseqrev = reverse($newseq);
        my $newseqrevcom = $newseqrev;
        $newseqrevcom =~tr/acgtnACGTN/tgcanTGCAN/;
	if ($desc) {
		if ($desc=~/^\s*1\b/) {
			$desc=~s/1\b/2/;
		} elsif ($desc=~/^\s*2\b/) {
			$desc=~s/2\b/1/;
		}
	}
        print   $id,(($desc) ? " $desc" : ''),"\n",
                $newseqrevcom,"\n",
                '+',"\n",
                $newqualrev,"\n";
    }
}
close(IN);

