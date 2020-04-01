#!/usr/bin/env perl

use strict;
use warnings;

use Storable;

use DBI;

my $data_source='dbi:mysql:database=pgpm;host=hammer.fcav.unesp.br';
my $username='dgpinheiro';
my $auth='b101nf0';

my $dbh = DBI->connect($data_source, $username, $auth, { RaiseError => 1, AutoCommit => 0 });


my $infile=$ARGV[0];

die "Missing input file (.tree)" unless ($infile);
die "Wrong input file ($infile)" unless (-e $infile);

my $annfile=$ARGV[1];

die "Missing annotation file (.txt)" unless ($annfile);
die "Wrong input file ($annfile)" unless (-e $annfile);


my $hr_annot;

if (-e "$annfile.dump") {
	$hr_annot = retrieve("$annfile.dump");
} else {
	open(ANN, "<", $annfile) or die $!;

	while(my $l = <ANN>) {
		chomp($l);
		my ($id, $desc)=split(/\t/, $l);
		$desc=~s/\sn=\d+ .*$//;
		#print $id,"\t",$desc,"\n";
		$hr_annot->{$id} = $desc;
	}
	close(ANN);

	store $hr_annot, "$annfile.dump";
}

my $outfile=$ARGV[2];
die "Missing output file (text file with annotation)" unless ($outfile);

my $annsf=$ARGV[3];

my %sf;
if ($annsf) {
	die "Wrong subfamily annotation file (last.subfam)" unless (-e $annsf);

	open(SFA, "<", $annsf) or die $!;
	my $sfid;
	while (<SFA>) {
		chomp;
		if ($_=~/\%subfamily (\S+)/) {
			$sfid=$1;
		} elsif (($sfid)&&($_=~/^>(\S+)/)) {
			my $acc=$1;
			$sf{ $acc } = $sfid;
		}
	}
	close(SFA);
}

my $seed;
my $seedfa=$ARGV[5];
if ($seedfa) {
	die "Wrong seed fasta file (seed.fa)" unless (-e $seedfa);

	open(SEED, "<", $seedfa) or die $!;
	while (<SEED>) {
		chomp;
		if ($_=~/^>(\S+)/) {
			$seed = $1;
		}
	}
	close(SEED);

	die "Not found any seed in seed fasta file ($seedfa)" unless ($seed);
}

my %msa;
my $msafa=$ARGV[4];
if ($msafa) {
	die "Wrong multiple alignment fasta file (msa.fa)" unless (-e $msafa);

	open(MSA, "<", $msafa) or die $!;
	while (<MSA>) {
		chomp;
		if ($_=~/^>(\S+)/) {
			my $sid=$1;
			$msa{$sid} = undef;
		}
	}
	close(MSA);
}

open(IN, "<", $infile) or die $!;
open(OUT, ">", $outfile) or die $!;

print OUT "#taxa\tdescription\tsubfamily\n" if ($annsf);

my $hr_tax;
my %tax;

while(my $l = <IN>) {
	while($l=~/\b([A-Z][A-Za-z_\.0-9]+)\b/g) {
		my $id = $1;
		my $uniref=$id;
		$uniref=~s/UniRef100_//;
		$hr_annot->{$id}||='';
		print OUT $id,"\t","$uniref $hr_annot->{$id}";
		if (exists $msa{ $id }) {
			print OUT '(','*';
			if ($seed eq $id) {
				print OUT '*';
			}
			print OUT ')';
		}
		
		if (exists $hr_tax->{ $id }) {
			print OUT ' / '.$tax{ $hr_tax->{ $id } };
		} else {
			my ($common_tax_id) = $dbh->selectrow_array("SELECT common_tax_id FROM uniref WHERE uniref_id=?", undef, $id);
			unless ($common_tax_id) {
				die "Not found common_tax_id FOR uniref_id=$id";
			}
			$hr_tax->{ $id } = $common_tax_id;
			unless (exists $tax{ $common_tax_id }) {
				$tax{ $common_tax_id } = `GetTaxInfo.pl $common_tax_id | tail -1 | cut -f 5`;
			}
		}
		
		if ($annsf) {
			print OUT " [$sf{$id}]\t",$sf{ $id };
		} else {
			print OUT "\t",'';
		}
		print OUT "\n";
	}
}

close(IN);
close(OUT);

$hr_annot = undef;
