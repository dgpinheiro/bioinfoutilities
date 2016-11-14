#!/bin/bash

infile=$1
outfile=$2

if [ ! ${infile} ]; then
	echo "Missing input file"
	exit
else
	if [ ! -e ${infile} ]; then
		echo "Wrong input file (${infile})"
		exit
	fi
fi

infbn=`basename ${infile} .gff`;

if [ ! ${outfile} ]; then
	outfile="${infbn}_fixed.gff"
fi

grep -v -P '\t(cDNA_match|match|repeat_region|region|D_loop|sequence_feature|binding_site)\t' ${infile} | \
perl -F"\t" -lane '
	INIT { 
		my %gene; 
		my %rna;
	} 
	if ($_=~/^#/) { 
		print $_; 
	} else { 
		if ($F[2] eq "gene") { 
			my ($ID)=$F[8]=~/ID=([^;]+)/; 
			$gene{$ID}=[@F]; 
		} elsif ( ($F[2]=~/[tr]RNA/) && ($F[8]!~/Parent/) ) { 
			my $eight = $F[8]; 
			$F[8] =~ s/ID=([^;]+)/ID=gene_$1/; 
			my $id=$1; 
			my $two = $F[2]; 
			$F[2] = "gene"; 
			print join("\t", @F); 
			$F[2]="RNA"; 
			$F[8]=$eight; 
			$F[8].=";Parent=gene_$id"; 
		} elsif ($F[2]=~/\b(?:(?:t?m|nc|SRP_|lnc_|t|r)RNA|primary_transcript|transcript|\w?_gene_segment)\b/) { 
			$F[2]="RNA"; 
			my ($Parent)=$F[8]=~/Parent=([^;]+)/; 
			$gene{$Parent}=undef;
		} elsif ($F[2]=~/\b(CDS|exon)\b/){  
			my ($Parent)=$F[8]=~/Parent=([^;]+)/; 
			if (exists $gene{$Parent}) { 
				my ($gene_id) = $Parent; die "Not found id for gene $Parent" unless ($gene_id);
				$gene{$gene_id} = undef;

				unless (exists $rna{ "rna_$gene_id" }) { $rna{ "rna_$gene_id" } = [@F]; $rna{ "rna_$gene_id" }->[2] = "RNA"; $rna{ "rna_$gene_id" }->[8] =~s/ID=[^;]+/ID=rna_$gene_id/; } else {
					if ($rna{ "rna_$gene_id" }->[4] < $F[4]) {
						$rna{ "rna_$gene_id" }->[4] = $F[4];
					}
					if ($rna{ "rna_$gene_id" }->[3] > $F[3]) {
						$rna{ "rna_$gene_id" }->[3] = $F[3];
					}
				}

				$F[8]=~s/Parent=([^;]+)/Parent=rna_$gene_id/; 

			}  
		} elsif ($F[2]=~/\b(miRNA)\b/){ 
			$F[2]=$1; 
		}  
		print join("\t", @F); 
	} 
	END {  
		foreach my $k (keys %rna) { 
			print join("\t", @{ $rna{$k} } ); 
		} 
		foreach my $k (keys %gene) { 
			if (defined $gene{$k}) { 
				$gene{$k}->[2]="RNA"; 
				$gene{$k}->[8]=~s/ID=([^;]+)/ID=rna_$1/; 
				$gene{$k}->[8].=";Parent=$k"; 
				print join("\t", @{$gene{$k}}); 
				$gene{$k}->[2]="exon"; 
				$gene{$k}->[8]=~s/ID=([^;]+)/ID=exon_$1/; 
				$gene{$k}->[8]=~s/Parent=([^;]+)/Parent=rna_$1/; 
				print join("\t", @{$gene{$k}}); 
			} 
		} 
	} ' | \
fixGFFdupID.pl | \
sortGFF.pl > ${outfile}

