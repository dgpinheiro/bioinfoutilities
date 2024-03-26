#!/bin/bash

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1 
#Diretorio que contem os arquivos processados

reftrinity=$2
# Arquivo da montagem de novo Trinity.fa

refprot=$3
# Diretorio que contem o proteoma de referencia mais proximo para ser utilizado como base

# output - diretório para armazenar o resultado do processo de montagem
output=$4

if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input directory ${input}"
                exit
        fi
fi

if [ ! ${reftrinity} ]
then
        echo "Missing assembled fasta file"
        exit
else
        if [ ! -e ${reftrinity} ]
        then
                echo "Wrong assembled fasta file  ${reftrinity}"
                exit
        fi
fi

if [ ! ${refprot} ]
then
        echo "Missing proteomic reference fasta file"
        exit
else
        if [ ! -e ${refprot} ]
        then
                echo "Wrong proteomic reference fasta file ${refprot}"
                exit
        fi
fi


if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi


num_threads=4
 

mkdir -p ${output}/Transrate_good

basedir_out="${output}/Transrate_good"


leftreads=()
rightreads=()

leftreads=( $( find ${input} -name '*_1.fastq' ) )
for r in ${leftreads[@]}; do
	r=`echo ${r} | sed 's/_1\./_2\./'`
	rightreads=(${rightreads[@]} ${r})
done

echo -e "Runnig Transrate"

transrate --reference ${refprot} \
          --assembly ${reftrinity} \
          --left $(IFS=, ; echo "${leftreads[*]}") \
          --right $(IFS=, ; echo "${rightreads[*]}") \
          --threads ${num_threads} \
          --output ${basedir_out} \
	  > ${output}/Transrate.log.out.txt \
	  2> ${output}/Transrate.log.err.txt
