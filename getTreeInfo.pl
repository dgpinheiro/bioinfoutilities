#!/usr/bin/env perl

use strict;
use warnings;

use Storable;

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
        my ($taxid) = $desc=~/TaxID=(\d+)/;
		$desc=~s/\sn=\d+ .*$//;
		#print $id,"\t",$desc,"\n";
		$hr_annot->{$id}->{'desc'} = $desc;
		$hr_annot->{$id}->{'taxid'} = $taxid;
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

print OUT "#taxa\tdescription\tsubfamily\ttaxonomy\ttax_id\n" if ($annsf);

my $hr_tax;
my %tax;

while(my $l = <IN>) {
	while($l=~/\b([A-Z][A-Za-z_\.0-9]+)\b/g) {
		my $id = $1;
		my $uniref=$id;
		$uniref=~s/UniRef100_//;
		$hr_annot->{$id}->{'desc'}||='';
		print OUT $id,"\t","$uniref $hr_annot->{$id}->{'desc'}";
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
			my ($common_tax_id);
            my $try=1;
            while ( (! -e "/tmp/$id.tab") && ($try<=3 ) ) {
                `curl -s https://www.uniprot.org/uniref/$id.tab > /tmp/$id.tab`;
                $try++;
                sleep(3);
            }

            open(TST, "<", "/tmp/$id.tab") or die "File not found /tmp/$id.tab";
            my $first=<TST>;
            chomp($first);
            close(TST);

            if ($first=~/^<!DOCTYPE/) {
                $common_tax_id=$hr_annot->{$id}->{'taxid'} if ($hr_annot->{$id}->{'taxid'});
            } else {
                $common_tax_id = `cat /tmp/$id.tab | cut -f 6  | tail -n+2 | GetLCA.pl`;
                chomp($common_tax_id);
                $common_tax_id||=1;
            }
            
			unless ($common_tax_id) {
				die "Not found common_tax_id FOR uniref_id=$id";
			}
            if ($common_tax_id !~ /^\d+$/) {
                die "This is not a number $common_tax_id";
            }
            $common_tax_id+=0;
            
			$hr_tax->{ $id } = $common_tax_id;
			unless (exists $tax{ $common_tax_id }) {
                my $gettaxinfo_cmd = "GetTaxInfo.pl $common_tax_id | tail -1 | cut -f 5";
    		    $tax{ $common_tax_id } = `$gettaxinfo_cmd`;
                chomp($tax{ $common_tax_id });
			}
		}
		
		if ($annsf) {
            my $taxon = ( (($hr_tax->{$id})&&($tax{ $hr_tax->{$id} })) ? $tax{ $hr_tax->{$id} } : '');
            chomp($taxon);
			print OUT " [$sf{$id}] ".$taxon."\t",$sf{ $id },"\t",$taxon,"\t",$hr_tax->{$id};
		} else {
			print OUT "\t",'',"\t",'';
		}
		print OUT "\n";
	}
}

close(IN);
close(OUT);

$hr_annot = undef;
