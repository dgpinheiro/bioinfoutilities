#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

my $infile = $ARGV[0];
die "Missing infile file" unless ($infile);
die "Wrong infile file ($infile)" unless (-e $infile);

my $reffile = $ARGV[1];
die "Missing reference file(s)" unless ($reffile);

my %refdesc;
foreach my $rfile (glob("$reffile")) {
	print STDERR "Loading $rfile ...\n";
	
	if ($rfile =~ /\.gz$/) {
		open(REF, "-|", "zcat $rfile") or die $!;
	} else {
		open(REF, "-|", "cat $rfile") or die $!;
	}
	
	my $seqref = Bio::SeqIO->new(-fh=>\*REF, -format=>'FASTQ');
	
	while(my $seq=$seqref->next_seq()) {
		#print $seq->display_id(),"\t",$seq->desc(),"\n";
		$refdesc{ $seq->display_id() }->{ $seq->seq() } = $seq->desc();
	}
	
	close(REF);
}

if ($infile =~ /\.gz/) {
	open(IN, "-|", "zcat $infile") or die $!;
} else {
	open(IN, "-|", "cat $infile") or die $!;
}

my $seqin = Bio::SeqIO->new(-fh=>\*IN, -format=>'FASTQ');
my $seqout = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'FASTQ');

while(my $seq=$seqin->next_seq()) {
	#print $seq->display_id(),"\t",$refdesc{ $seq->display_id() }->{ $seq->seq() },"\n";
	if ((exists $refdesc{ $seq->display_id() }) && (exists $refdesc{ $seq->display_id() }->{ $seq->seq() } )) {
		$seq->desc( $refdesc{ $seq->display_id() }->{ $seq->seq() } );
	}
	$seqout->write_seq( $seq );
}

close(IN);
