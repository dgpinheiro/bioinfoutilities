#!/bin/bash

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

# validação do parâmetro "input"
if [ ! ${input} ]
then   
        echo "[ERROR] Missing input directory." 1>&2
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "[ERROR] Wrong input directory (${input})." 1>&2
                exit
        fi
fi

# output - diretório para armazenar o resultado do processo de montagem
output=$2

# validação do parâmetro "output"
if [ ! ${output} ]
then   
        echo "[ERROR] Missing output directory." 1>&2
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "[ERROR] Wrong output directory (${output})." 1>&2
                exit
        fi
fi

# Número de CORES para o processamento
# ATENÇÃO: Não exceder o limite da máquina
# Quantidade de memória para o processamento com Jellyfish
# ATENÇÃO: Não exceder o limite da máquina
THREADS=$3

if [ ! ${THREADS} ]; then
	THREADS=1
	echo "[WARNING] Using default value (${THREADS}) for the number of threads." 1>&2
	MEM=$3
else
	MEM=$4
fi


if [ ! ${MEM} ]; then
	MEM=10
	echo "[WARNING] Using default value (${MEM}) for the memory in Gb." 1>&2
fi

MEM="${MEM}G"

#echo -e "THREADS: ${THREADS}\nMEMORY: ${MEM}"

###
# Arquivos e diretórios de saída (output) 
#

basedir_out="${output}/"

renamed_out="${basedir_out}/renamed"

trinity_out="${basedir_out}/trinity_assembled"

mkdir -p ${renamed_out}
mkdir -p ${trinity_out}

left=()
left_singleton=()

right=()
right_singleton=()

echo "Performing renaming step ..."

for fastq in `ls ${input}/*.fastq`; do
	# obtendo nome do arquivo 
	fastqbn=`basename ${fastq}`;
	if [[ ! $fastqbn =~ \.bad_ ]]; then
		renamed_fastq="${renamed_out}/${fastqbn}"
		if [ ! -e ${renamed_fastq} ]; then
			echo -e "\tRenaming ${fastqbn} ..."
			if [[ ${fastqbn} =~ _1[\._] ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/1$/) { print $1"/1" } else { print $1 } } else if (NR%4==3) { print "+" } else { print $1 } }' ${fastq} > ${renamed_fastq}
			elif [[ ${fastqbn} =~ _2[\._]  ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/2$/) { print $1"/2" } else { print $1 } } else if (NR%4==3) { print "+" } else { print $1 } }' ${fastq} > ${renamed_fastq}
			else 
				echo "Warning: ${fastqbn} discarded!"
			fi
		fi
		
		if [[ ${fastqbn} =~ _1[\._] ]]; then
			if [[ ${fastqbn} =~ singletons ]]; then
				if [ -s ${renamed_fastq} ]; then
					left_singleton=($(printf "%s\n" ${left_singleton[@]} ${renamed_fastq} | sort -u ))
				fi
			else
				left=($(printf "%s\n" ${left[@]} ${renamed_fastq}  | sort -u ))
			fi
		elif [[ ${fastqbn} =~ _2[\._] ]]; then
			if [[ ${fastqbn} =~ singleton ]]; then
				if [ -s ${renamed_fastq} ]; then
					right_singleton=($(printf "%s\n" ${right_singleton[@]} ${renamed_fastq}  | sort -u ))
				fi
			else
				right=($(printf "%s\n" ${right[@]} ${renamed_fastq}  | sort -u ))
			fi
		else
			echo "Warning: ${fastqbn} discarded!"
		fi
	fi
done


#for l in ${left[@]}; do
#	echo -e "L: ${l}";
#done
#
#for r in ${right[@]}; do
#	echo -e "R: ${r}";
#done
#
#for ls in ${left_singleton[@]}; do
#	echo -e "LS: ${ls}";
#done
#
#for rs in ${right_singleton[@]}; do
#	echo -e "RS: ${rs}";
#done

if [ ! -e ${trinity_out}/Trinity.fasta ]; then
	
	echo -e "Assembling step (Trinity) ..."
	
	rm -fr ${trinity_out}
	mkdir -p ${trinity_out}
	
	Trinity --output ${trinity_out}/ \
		--seqType fq \
		--max_memory ${MEM} \
		--CPU ${THREADS} \
		--min_per_id_same_path 98 \
		--path_reinforcement_distance 25 \
		--group_pairs_distance 200 \
		--min_kmer_cov 5 \
		--min_glue 4 \
		--min_contig_length 300 \
		--SS_lib_type RF \
		--left $(IFS=, ; echo "${left[*]},${left_singleton[*]}") \
		--right $(IFS=, ; echo "${right[*]},${right_singleton[*]}") \
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt

	mv ${basedir_out}/trinity_assembled.Trinity.fasta ${trinity_out}/Trinity.fasta
	mv ${basedir_out}/trinity_assembled.Trinity.fasta.gene_trans_map ${trinity_out}/Trinity.fasta.gene_trans_map

fi
