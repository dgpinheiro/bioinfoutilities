#!/bin/bash

infile=$1
indir=$2
outdir=$3

if [ ! ${infile} ]; then
	echo "[ERROR] Missing input data file!" 2>&1
	exit
else
	if [ ! -e ${infile} ]; then
		echo "[ERROR] Wrong input data file (${infile})!" 2>&1
		exit
	fi
fi

if [ ! ${indir} ]; then
	echo "[ERROR] Missing input raw directory!" 2>&1
	exit
else
	if [ ! -d ${indir} ]; then
		echo "[ERROR] Wrong input raw directory (${indir})!" 2>&1
		exit
	fi
fi

if [ ! ${outdir} ]; then
	echo "[ERROR] Missing output symlink directory!" 2>&1
	exit
else
	if [ ! -d ${outdir} ]; then
		echo "[ERROR] Wrong output symlink directory (${outdir})!" 2>&1
		exit
	fi
fi

curdir=`pwd`
absindir=`readlink -f ${indir}`

while IFS= read -r line;
do
   oldfile=`echo ${line} | cut -d ' ' -f 1 `;
   newfile=`echo ${line} | cut -d ' ' -f 2 `;

   cd ${outdir}

   if [ ! -e ${absindir}/${oldfile} ]; then
	   echo "[ERROR] Not found ${absindir}/${oldfile}!" 2>&1
	   exit
   fi
   ln -s ${absindir}/${oldfile} ${newfile}
   
   cd ${curdir}

done < ${infile}
