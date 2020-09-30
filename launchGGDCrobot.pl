#!/usr/bin/perl

use strict;
use warnings;

# Limite estabelecido na ferramenta online
use constant REFLIMIT=>75;

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
    my $ds=0;
    my $c=1;
    for (my $g=$f+1; $g<=$#file; $g++) {
        if ($c > REFLIMIT) {
            $c=1;
            $ds++;
        }
        push(@{$other[$ds]}, $file[$g]);
        $c++;
    }
    foreach my $dsc (0..$#other) {
        #print $file[$f] . ' x ' . join(",", @{ $other[$dsc] })." (Dataset ".($dsc+1).")\n";
        print $file[$f] . ' x ' . scalar(@{ $other[$dsc] })." genomes (Dataset ".($dsc+1).")\n";
        my $cmd="./robot.pl -q $file[$f] -r ".join(",", @{ $other[$dsc] }). ' -e '.$email;
        print `$cmd`;
    }
}
