#!/bin/bash

abddir=$1
transcriptome=$2
outdir=$3
number=$4

if [ ! ${number} ]; then
	number=1000
	echo "[WARNING] Using default ${number} reads!" 2>&1
fi

declare -a abdlist

abdlist=()
if [ ! ${abddir} ]; then
	echo "[ERROR] Missing abundance input directory!" 2>&1
	exit
else
	if [ ! -d ${abddir} ]; then
		echo "[ERROR] Wrong abundance input directory (${abddir})!" 2>&1
		exit
	fi

	for abdfile in `ls ${abddir}/abundance_?.txt`; do
		abdlist+=(${abdfile})
	done

	if [ ${#abdlist[@]} == 0 ]; then
		echo "[ERROR] Not found abundance_?.txt files in ${abddir}" 2>&1
		exit
	fi
fi

if [ ! ${transcriptome} ]; then
	echo "[ERROR] Missing transcriptome file!" 2>&1
	exit
else
	if [ ! -e ${transcriptome} ]; then
		echo "[ERROR] Wrong transcriptome file (${transcriptome})!" 2>&1
		exit
	fi
fi

if [ ! ${outdir} ]; then
	echo "[ERROR] Missing output directory!" 2>&1
	exit
else
	if [ ! -d ${outdir} ]; then
		echo "[ERROR] Wrong output directory (${outdir})!" 2>&1
		exit
	fi
fi

for abdfile in ${abdlist[@]}; do
	samplen=`echo ${abdfile} | sed 's/^.*abundance_\([A-Z]\)\.txt/\1/'`;
	echo "Generaring simulation for group ${samplen} ..."
	for rep in `seq 1 3`; do
		echo -e "\tReplicate ${rep}"
		simLib.pl -a ${abdfile} -i ${transcriptome} -rf /usr/local/bioinfo/simNGS/data/s_4_0099.runfile -n ${number} -o ${outdir} -p SAMPLE${samplen}${rep}
	done
done
