infile=$1
prefix=$2
maturefile=$3
gffile=$4

if [ ! ${infile} ]; then
	echo "Missing input file"
	exit
else
	grp=`echo ${infile} | grep -P '\.fa$'`;
	if [ ! ${grp} ]; then
		echo "The input file needs to have .fa suffix"
		exit
	fi
fi

if [ ! ${prefix} ]; then
	echo "Missing prefix"
	exit
fi

if [ ! ${maturefile} ]; then
	echo "Missing mature file"
	exit
fi

if [ ! ${gffile} ]; then
	echo "Missing gff"
	exit
fi

bn=`basename ${infile} .fa`;
bngff=`basename \`basename ${gffile} .gff \` .gff3`

pita_prediction.pl -utr ${infile} -mir ${maturefile} -prefix ${prefix} -l 6-8 -gu '6;0,7;1,8;1' -m '6;0,7;0,8;1' > ${prefix}_pita_output.txt 2>&1

./selectPITA.pl -i ${prefix}_pita_results_targets.tab -t -10 > Targets_${bn}_s-10.tmp

./putGeneInfo.pl -i Targets_${bn}_s-10.tmp -g ${gffile} -c 1 -r '^[^\|]+\|[^\|]+\|([^\|]+)\|' > Targets_${bn}_s-10-${bngff}.txt

