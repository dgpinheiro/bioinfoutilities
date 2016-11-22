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

grep -v -P '\t(cDNA_match|match|repeat_region|region|D_loop|sequence_feature|binding_site|mobile_genetic_element)\t' ${infile} | \
fixGFFgenestruct.pl | \
fixGFFdupID.pl | \
sortGFF.pl > ${outfile}
