infa=$1
mirmature=$2
output=$3

if [ ! ${infa} ]; then
	echo "Missing input fasta file"
	exit
else
	if [ ! -f ${infa} ]; then
		echo "Wrong input fasta file (${infa})"
		exit
	fi	
fi

if [ ! ${mirmature} ]; then
	echo "Missing mature fasta file"
	exit
else
	if [ ! -f ${mirmature} ]; then
		echo "Wrong mir mature fasta file (${mirmature})"
		exit
	fi	
fi

if [ ! ${output} ]; then
	echo "Missing output dir"
	exit
else 
	if [ ! -d ${output} ]; then
		echo "Wrong output dir (${output})"
		exit
	fi	
fi

mkdir -p ${output}/calibration
mkdir -p ${output}/targets

inbase=`basename $infa`

RNAcalibrate -f 2,7 -k 5000 -t ${infa} -q ${mirmature} -m 100000 > ${output}/calibration/calibrate_${inbase}.txt

perl runRNAhybrid.pl -q ${mirmature} -t ${infa} -o ${output}/targets/ -r OUT_ -p ${output}/calibration/calibrate_${inbase}.txt -c ' -f 2,7 -c -e -20 -p 0.1 -g png -m 100000'

