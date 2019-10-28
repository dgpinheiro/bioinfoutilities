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

grep -v -P '\t(cDNA_match|match|inverted_repeat|repeat_region|region|D_loop|sequence_feature|binding_site|mobile_genetic_element|origin_of_replication|STS|sequence_alteration|Shine_Dalgarno_sequence|ribosome_entry_site|recombination_feature|stem_loop|telomere|centromere|riboswitch|minus_35_signal|minus_10_signal|protein_binding_site|golden_path_region|match_part|transcription_start_site|match|exon_junction|TF_binding_site|BAC_cloned_genomic_insert|DNA_motif|RNAi_reagent|TF_binding_site|breakpoint|chromosome|chromosome_band|complex_substitution|deletion|enhancer|exon_junction|five_prime_UTR|golden_path_region|insertion_site|insulator|intron|match|match_part|modified_RNA_base_feature|oligo|oligonucleotide|origin_of_replication|orthologous_region|orthologous_to|pcr_product|point_mutation|polyA_site|protein|protein_binding_site|region|regulatory_region|repeat_region|rescue_fragment|sequence_variant|silencer|syntenic_region|tandem_repeat|three_prime_UTR|transcription_start_site|transposable_element|transposable_element_insertion_site|uncharacterized_change_in_nucleotide_sequence|direct_repeat)\t' ${infile} | 
grep -v -P '^[ACGTN]+$' | 
grep -v '^>' |
grep -v -P '^(##sequence-region|##species|#!genome-build|!processor|#!gff-spec-version)' | \
fixGFFgenestruct.pl | \
fixGFFdupID.pl | \
sortGFF.pl > ${outfile}
