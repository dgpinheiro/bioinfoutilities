#!/usr/bin/env perl
 
use strict;
use warnings;
use Data::Dumper;
 
use File::Basename;
use Bio::FeatureIO;
use FileHandle;
use POSIX 'isatty';

my $infile=$ARGV[0];

my $fh;

if ($infile) {

    die "Wrong input file ($infile)" unless (-e $infile);

    $fh = FileHandle->new;
    $fh->open("<$infile");

} else {
    unless (isatty(*STDIN)) {
        $fh = \*STDIN;
    } else {
        die "Missing input file or STDIN data";
    }
}

my %parent;

my @data;

my %GeneID;

while(<$fh>) {
	if ($_=~/^#/) {
		next;
	} else {
		chomp;
		my (@F) = split(/\t/, $_);
		
		if ($F[2] =~ /exon|CDS/i) {
			push(@data, \@F);
		} else {
			my ($ID)=$F[8]=~/ID=([^;]+)\b/;
			die "Missing ID ($_)" unless ($ID);
			if ($F[2] ne 'gene') {
				my ($Parent)=$F[8]=~/Parent=([^;]+)\b/;
				die "Missing Parent ($_)" unless ($Parent);
				$parent{ $ID } = $Parent;
			} else {
				my ($GID)=$F[8]=~/GeneID:([^;,]+)\b/;
				unless ($GID) {
					my ($gene)=$F[8]=~/gene=([^;]+)\b/;	
					if ($gene) {
						($GID)=$gene=~/LOC(\d+)/;
					}
				}
				if ($GID) {
					$GeneID{ $ID } = $GID;
				}
			}
		}
	}
}

$fh->close();
	
foreach my $ar_data (@data) {
	my @F = @{ $ar_data };
	my ($ID)=$F[8]=~/ID=([^;]+)\b/;
	die "Missing ID ($_)" unless ($ID);
	my ($Parent)=$F[8]=~/Parent=([^;]+)\b/;
	die "Missing Parent ($_)" unless ($Parent);
	my ($gene)=$F[8]=~/gene=([^;]+)\b/;	
	
	my ($transcript_id, $gene_id);
	$transcript_id=$Parent;
	$gene_id=((exists $GeneID{ $parent{ $transcript_id } }) ? $GeneID{ $parent{ $transcript_id } } : $parent{ $transcript_id }); 
	
	$F[8]="gene_id \"$gene_id\"; ".(($gene) ? "gene_name=\"$gene\"; " : "")."transcript_id \"$transcript_id\";";
	print join("\t", @F),"\n";
}

