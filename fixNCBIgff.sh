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

grep -v -P '\t(cDNA_match|match|repeat_region|region|D_loop|sequence_feature|binding_site|mobile_genetic_element|origin_of_replication|STS|sequence_alteration|Shine_Dalgarno_sequence|ribosome_entry_site|recombination_feature|stem_loop|RNase_P_RNA|telomere|centromere|telomerase_RNA|RNase_MRP_RNA|riboswitch|minus_35_signal|minus_10_signal|protein_binding_site)\t' ${infile} | 
grep -v -P '^(##sequence-region|##species|#!genome-build|!processor|#!gff-spec-version)' | \
fixGFFgenestruct.pl | \
fixGFFdupID.pl | \
sortGFF.pl > ${outfile}
