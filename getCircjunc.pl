#!/usr/bin/env perl

use strict;
use warnings;

my $infile=$ARGV[0];

my $size=$ARGV[1];

$size||=100;

die "Missing input file" unless ($infile);
die "Wrong input file ($infile)" unless (-e $infile);

use Bio::SeqIO;

my $seqin = Bio::SeqIO->new(-file=>$infile, -format=>'FASTA');

my $seqout = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'FASTA');

while(my $seq = $seqin->next_seq()) {

    my $seqobj = Bio::PrimarySeq->new (
              -seq              => $seq->subseq( $seq->length()-$size+1, $seq->length()).$seq->subseq(1,$size),
              -id               => $seq->display_id().'_circ_junc',
              -alphabet         => 'dna'
              );
    print STDERR "Get circular junction size: ".$seqobj->length(),"\n";
    $seqout->write_seq($seqobj);
}
