#!/bin/bash

# input - diretório contendo os arquivos de entrada no formato .fastq (renomeados)
input=$1

# validação do parâmetro "input"
if [ ! ${input} ]
then   
        echo "Missing input (renamed for Trinity) directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input (renamed for Trinity) directory ${input}"
                exit
        fi
fi

# trinity_output - diretório onde foi armazenado o resultado do processo de montagem
trinity_output=$2

# validação do parâmetro "trinity_output"
if [ ! ${trinity_output} ]
then   
        echo "Missing Trinity output directory"
        exit
else   
        if [ ! -d ${trinity_output} ]
        then   
                echo "Wrong Trinity output directory ${trinity_output}"
                exit
        fi
fi

# output - diretório para armazenar o resultado da avaliação de abundância
output=$3

# validação do parâmetro "output"
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

num_threads="4"

###
# Arquivos e diretórios de saída (output) 
#

abundance_out="${output}/abundance"
mkdir -p ${abundance_out}

left=()
right=()

echo "Collecting reads step ..."


left=($(find ${input} -type f -name '*_1.fastq'))

rm -f ${abundance_out}/samples.txt
rm -f ${abundance_out}/quant_files.txt
rm -f ${abundance_out}/groups.txt

#echo -e "id\tname\tgroup" > ${abundance_out}/groups.txt

for l in ${left[@]}; do
	repname=`basename ${l} | sed 's/\..*$//'`
	condname=`echo ${repname} | sed 's/_B[0-9]\+//'`
	r=`echo ${l} | sed 's/_1.fastq/_2.fastq/'`
	right=(${right[@]} ${r})

	echo -e "${condname}\t${abundance_out}/${repname}\t${l}\t${r}" >> ${abundance_out}/samples.txt
	# Se for utilizado kallisto é necessário alterar no nome quant.sf para abundances.tsv
	echo -e "${abundance_out}/${repname}/quant.sf" >> ${abundance_out}/quant_files.txt

	#echo -e "${repname}\t${repname}\t${condname}" >> ${abundance_out}/groups.txt
	echo -e "${condname}\t${repname}" >> ${abundance_out}/samples_described.txt
done


echo "Left files..: ${left[*]}"
echo "Right files.: ${right[*]}"

trinity_fasta=`find ${trinity_output} -type f -name 'Trinity.fasta'`
trinity_trans_map=`find ${trinity_output} -type f -name Trinity.fasta.gene_trans_map`

echo "Using transcriptome assembly: ${trinity_fasta} and ${trinity_trans_map}"

echo "Estimating abundances ..."

							#--est_method	kallisto \
${TRINITY_HOME}/util/align_and_estimate_abundance.pl 	--transcripts	${trinity_fasta} \
							--samples_file	${abundance_out}/samples.txt \
							--est_method	salmon \
							--gene_trans_map ${trinity_trans_map} \
							--prep_reference \
							--thread_count ${num_threads} \
							--seqType fq \
							--SS_lib_type RF \
							--output_dir ${abundance_out} \
							 > ${abundance_out}/align_and_estimate_abundance.log.out.txt \
							2> ${abundance_out}/align_and_estimate_abundance.log.err.txt


echo "Constructing abundance matrix ..."

							#--est_method	kallisto \
${TRINITY_HOME}/util/abundance_estimates_to_matrix.pl	--est_method salmon \
							--gene_trans_map ${trinity_trans_map} \
							--name_sample_by_basedir \
							--cross_sample_norm none \
							--quant_files ${abundance_out}/quant_files.txt \
							--out_prefix ${abundance_out}/abundance \
							 > ${abundance_out}/abundance_estimates_to_matrix.log.out.txt \
							2> ${abundance_out}/abundance_estimates_to_matrix.log.err.txt 


echo "Calculating Differentially Expressed Genes ..."

mkdir -p ${abundance_out}/gene/DESeq2

${TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl \
          				    --matrix ${abundance_out}/abundance.gene.counts.matrix \
          				    --method DESeq2 \
          				    --samples_file ${abundance_out}/samples_described.txt \
					    --output ${abundance_out}/gene/DESeq2 \
					    --reference_sample TREAT
mkdir -p ${abundance_out}/gene/edgeR

${TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl \
          				    --matrix ${abundance_out}/abundance.gene.counts.matrix \
          				    --method edgeR \
          				    --samples_file ${abundance_out}/samples_described.txt \
					    --output ${abundance_out}/gene/edgeR \
					    --reference_sample TREAT

#run-DESeq2.R 	--in="${abundance_out}/abundance.gene.counts.matrix"  \
#		--groups="${abundance_out}/groups.txt" \
#		--out="${abundance_out}/DEG" \
#		 > ${abundance_out}/DEG/run-DESeq2.log.out.txt \
#		2> ${abundance_out}/DEG/run-DESeq2.log.err.txt


