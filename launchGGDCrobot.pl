#!/usr/bin/perl

use strict;
use warnings;

my $dataf = $ARGV[0];
my $email = $ARGV[1];

die "Missing file with list of genome fasta file paths" unless ($dataf);
die "Wrong file with list of genome fasta file paths ($dataf)" unless (-e $dataf);

$email||='bioinfo.fcav@gmail.com';

print STDERR "The results will be send to $email.\n";

my @file;
open(IN, "<", $dataf) or die $!;
while(<IN>) {
    chomp;
    $_=~s/^\s+//;
    $_=~s/\s+$//;
    push(@file, $_);
}
close(IN);

for (my $f=0; $f<$#file; $f++) {
    my @other;
    for (my $g=$f+1; $g<=$#file; $g++) {
        push(@other, $file[$g]);
    }
    #print $file[$f] . ' x ' . join(",", @other)."\n";
    my $cmd="GGDCrobot.pl -q $file[$f] -r ".join(",", @other). ' -e '.$email.' ';
    print `$cmd`;
}
