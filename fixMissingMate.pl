#!/usr/bin/perl

use strict;
use warnings;

#use File::Temp qw/ tempfile tempdir /;

my $fastqin=$ARGV[0];

die "Missing input file (.fastq)" unless ($fastqin);

my $left_multimapped=$ARGV[1];
die "Missing left_multimapped (.bam)" unless ($left_multimapped);
die "Wrong left_multimapped (.bam) ($left_multimapped)" unless (-e $left_multimapped);

my $right_multimapped=$ARGV[2];
die "Missing right_multimapped (.bam)" unless ($right_multimapped);
die "Wrong right_multimapped (.bam) ($right_multimapped)" unless (-e $right_multimapped);

#my $dir = tempdir( CLEANUP => 1 );
#my ($fh, $filename) = tempfile( DIR => $dir );

my @seq;
my $c = 0;
foreach my $fqin (split(/,/, $fastqin)) {

    die "Wrong input file (.fastq) ($fqin)" unless (-e $fqin);

    open(IDX, $fqin) or die $!;
    while(<IDX>) {
        chomp;
        if ($.%4 == 1) {
            my ($name) = $_;
            $name=~s/^@//;
            $name=~s/ .*$//;
            push(@seq, $name);
        }            
    }
    close(IDX);
}

open(BAM, "-|", "samtools view -H $left_multimapped") or die $!;
while(<BAM>) {
    print $_;
}
close(BAM);

my %bam_file = (   'left'=>{'file'=>$left_multimapped, 'bit'=>6},
                   'right'=>{'file'=>$right_multimapped, 'bit'=>7} );

foreach my $k ( keys %bam_file ) {

    my $multimapped_file = $bam_file{ $k }->{'file'};
    my $bit = $bam_file{ $k }->{'bit'};

    open(BAM, "-|","samtools view $multimapped_file") or die $!;
    while(<BAM>) {
        chomp;
        my @F=split(/\t/, $_);
        die "Not found sequence ID $F[0]" unless ($seq[$F[0]-1]);
        $F[0] = $seq[$F[0]-1];
        
        my $revflag=reverse(sprintf("%012b", $F[1]));
        substr($revflag, 0, 1, 1);
        substr($revflag, 3, 1, 1);
        substr($revflag, $bit, 1, 1);
        substr($revflag, 8, 1, 0);
        my $adjflag=reverse($revflag);
        $F[1] = oct('0b'.$adjflag);
        
        print join("\t", @F),"\n";
    }
    close(BAM);
}

