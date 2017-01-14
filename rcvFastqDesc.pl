#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

my $infile = $ARGV[0];
die "Missing infile file" unless ($infile);
die "Wrong infile file ($infile)" unless (-e $infile);

my $reffile = $ARGV[1];
die "Missing reffile file" unless ($reffile);
die "Wrong reffile file ($reffile)" unless (-e $reffile);

my $seqout = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'FASTQ');
my $seqref;
if ($reffile =~ /\.gz$/) {
	$seqref = Bio::SeqIO->new(-file=>"gunzip -c $reffile |", -format=>'FASTQ');
} else {
	$seqref = Bio::SeqIO->new(-file=>"$reffile", -format=>'FASTQ');
}
my %refdesc;
while(my $seq=$seqref->next_seq()) {
	#print $seq->display_id(),"\t",$seq->desc(),"\n";
	$refdesc{ $seq->display_id() }->{ $seq->seq() } = $seq->desc();
}
my $seqin;
if ($infile =~ /\.gz/) {
	$seqin = Bio::SeqIO->new(-file=>"gunzip -c $infile |", -format=>'FASTQ');
} else {
	$seqin = Bio::SeqIO->new(-file=>$infile, -format=>'FASTQ');
}
my %ref;
while(my $seq=$seqin->next_seq()) {
	#print $seq->display_id(),"\t",$refdesc{ $seq->display_id() }->{ $seq->seq() },"\n";
	if ((exists $refdesc{ $seq->display_id() }) && (exists $refdesc{ $seq->display_id() }->{ $seq->seq() } )) {
		$seq->desc( $refdesc{ $seq->display_id() }->{ $seq->seq() } );
	}
	$seqout->write_seq( $seq );
}

