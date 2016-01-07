#!/usr/bin/perl

use strict;
use warnings;

my $fastqfile = $ARGV[0];

die "Missing single-end .fastq file" unless ($fastqfile);
die "Wrong single-end .fastq file ($fastqfile)" unless (-e $fastqfile);

my $id;
my $seq;
my $qual;

open(IN, "<", $fastqfile) or die $!;
while(<IN>) {
    chomp;
    if ($. % 4 == 1) {
        ($id)=$_=~/^(\S+)/;
    } elsif ($. % 4 == 2) {
        $seq=$_;
    } elsif ($. % 4 == 0) {
        $qual=$_;

        my $newqual=$qual;
        my $newseq=$seq;
        
        print   $id,"\n",
                $newseq,"\n",
                '+',"\n",
                $newqual,"\n";

        my $newqualrev = reverse($newqual);
        my $newseqrev = reverse($newseq);
        my $newseqrevcom = $newseqrev;
        $newseqrevcom =~tr/acgtnACGTN/tgcanTGCAN/;

        print   $id,"\n",
                $newseqrevcom,"\n",
                '+',"\n",
                $newqualrev,"\n";
    }
}
close(IN);

