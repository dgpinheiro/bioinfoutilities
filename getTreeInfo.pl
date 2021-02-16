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


$|=1;

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
            my ($id2)=$id=~/UniRef100_(\S+)/;

            $common_tax_id=&getCommonTaxId($hr_annot, $id, "https://www.uniprot.org/uniref/#ID#.tab", 6);
            $common_tax_id=&getCommonTaxId($hr_annot, $id2, "https://pir3.uniprot.org/uniref/?query=#ID#&fil=identity%3A1.0&columns=id%2Creviewed%2Corganism-id%2Cname&sort=score&format=tab", 3) unless ($common_tax_id);
            $common_tax_id=&getCommonTaxId($hr_annot, $id2, "https://www.uniprot.org/uniparc/?query=#ID#&columns=id%2Corganism-id&format=tab", 2) unless ($common_tax_id);

			unless ($common_tax_id) {
				die "Not found common_tax_id FOR uniref_id=$id";
			}

            $common_tax_id||=1;
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
			

sub getCommonTaxId {
    my ($hr_annot, $lid, $url, $ncut) = @_;
    my $try=1;
    $url=~s/#ID#/$lid/;
    while ( (! -e "/tmp/$lid.tab") && ($try<=3 ) ) {
        `(curl -sL "$url" > /tmp/$lid.tab) && (sync)`;
        $try++;
        sleep(3);
        `sync`;
    }
            
    my $first="";
    if ((-e "/tmp/$lid.tab")&&(! -z "/tmp/$lid.tab")) {
        open(TST, "<", "/tmp/$lid.tab") or die "File not found /tmp/$lid.tab";
        $first=<TST>;
        chomp($first);
        close(TST);
    } else {
        warn("Empty or non existent file: /tmp/$lid.tab");
        unlink("/tmp/$lid.tab");
        return(undef);
    }

    if ($first=~/^<!DOCTYPE/) {
        if ($hr_annot->{$lid}->{'taxid'}) {
            return($hr_annot->{$lid}->{'taxid'});
        } else {
            return(undef);
        }
    } else {
        my $ctid = `cat /tmp/$lid.tab | cut -f $ncut  | tail -n+2 | sed 's/;//g' | replace_merged_taxid.pl | GetLCA.pl`;
        chomp($ctid);
        if ($ctid) {
            if ($ctid !~ /^\d+$/) {
                warn "Not a number found for $lid";
                return(undef);
            }
        }            
        return($ctid);
    }
}
