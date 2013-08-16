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

mkdir -p ./calibration/

RNAcalibrate -f 2,7 -k 5000 -t ${infile} -q ${maturefile} -m 100000 > ./calibration/calibrate_${bn}.txt 2> ./calibration/calibrate_${bn}_err.txt

mkdir -p ./targets/

./runRNAhybrid.pl -q ${maturefile} -t ${infile} -o ./targets/ -r ${prefix} -p ./calibration/calibrate_${bn}.txt -c ' -f 2,7 -c -e -20 -p 0.1 -m 100000 '

./parseRNAhybridTargets.pl -i ./targets -e -20 -p 0.005 -g -s 2,7 > Targets_${bn}_e-20_p0.005_noseedGU2-7.tmp

./putGeneInfo.pl -i Targets_${bn}_e-20_p0.005_noseedGU2-7.tmp -g ${gffile} -c 1 -r '^[^\|]+\|[^\|]+\|([^\|]+)\|' > Targets_${bn}_e-20_p0.005_noseedGU2-7-${bngff}.txt

