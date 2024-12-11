#!/bin/bash

infile=$1

if [ ! ${infile} ]; then
	echo "[ERROR] Missing input fasta file" 2>&1
	exit
fi

sed 's/^>\(\S\+\).*/>\1/' $infile
