#!/bin/bash
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2012  Universidade de Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade de Estadual Paulista "Júlio de Mesquita Filho"
#  Faculdade de Ciências Agrárias e Veterinárias
#  Laboratório de Bioinformática (FCAV-BioInfo)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com


#############################
# Caminhos dos arquivos e diretórios de entrada (input - diretório com os arquivos .fastq)
# 

input="$1"
if [ ! ${input} ]
then   
	echo "Missing input path to directory containing *.fastq files"
        exit
else   
        if [ ! -d ${input} ]
        then   
		echo "Wrong input path (${input}). Please, check the first argument."
                exit
        fi
fi

outdir="$2"
if [ ! ${outdir} ]
then   
	outdir="."
else
	if [ ! -d ${outdir} ]
	then
		echo "Output directory (${outdir}) doesn't exist. Please, choose another one."
		exit
	fi	
fi


# Número de threads
threads=14

#############################
# Poda de adaptadores e primers
#

# Diretório base onde estão os adaptadores
adapt_barc_primer_path="./refs/adapter"

# scythe_adapters_3p_R1_path - caminho para o arquivo contendo os adaptadores 3p para R1
# scythe_adapters_3p_R2_path - caminho para o arquivo contendo os adaptadores 3p para R2
# (de acordo com o protocolo/tecnologia de sequenciamento) para o programa scythe
scythe_adapters_3p_R1_path="${adapt_barc_primer_path}/truseq_illumina_3p_adapt_R1.fa"
scythe_adapters_3p_R2_path="${adapt_barc_primer_path}/truseq_illumina_3p_adapt_R2.fa"

# cutadapt_adapters_5p_R1_path - caminho para o arquivo contendo os adaptadores 5p para R1
# cutadapt_adapters_5p_R2_path - caminho para o arquivo contendo os adaptadores 5p para R2
# (de acordo com o protocolo/tecnologia de sequenciamento) para o programa scythe
cutadapt_adapters_5p_R1_path="${adapt_barc_primer_path}/truseq_illumina_5p_adapt_R1.fa"
cutadapt_adapters_5p_R2_path="${adapt_barc_primer_path}/truseq_illumina_5p_adapt_R2.fa"

#########################################
###
# Arquivos e diretórios de saída (output) 
#

# diretório a partir de onde o script será executado
base_out="${outdir}"

processed_out="${base_out}/processed"
# fastqc_pre_out - caminho para o diretório onde serão criados os gráficos de qualidade dos dados brutos
fastqc_pre_out="${processed_out}/fastqc/pre"

# scythe_out - caminho para o diretório de saída do scythe
scythe_out="${processed_out}/scythe"

# cutadapt_out - caminho para o diretório de saída do cutadapt
cutadapt_out="${processed_out}/cutadapt"
	
# prinseq_out - caminho para o diretório de saída do prinseq
prinseq_out="${processed_out}/prinseq"

# fastqc_pos_out - caminho para o diretório onde serão criados os gráficos de qualidade dos dados processados
fastqc_pos_out="${processed_out}/fastqc/pos"


########################################################

echo "* Creating directories into base directory (${base_out})"

# Criando diretórios para as saídas dos programas que serão utilizados a seguir
mkdir -p ${fastqc_pre_out}
mkdir -p ${scythe_out}
mkdir -p ${cutadapt_out}
mkdir -p ${prinseq_out}
mkdir -p ${fastqc_pos_out}

###################################################
#FastQC
## percorrendo múltiplos valores retornados por um comando no shell (o comando irá listar todos os arquivos com extensão *_R1*.fastq

for fastq1 in `ls ${input}/*_R1*.fastq`; do
	fastq2=`echo ${fastq1} | sed 's/R1/R2/'`;
	
	# fqname1 - string com o nome do arquivo fastq
	fqname1=`basename ${fastq1} .fastq`
	# fqname2 - string com o nome do arquivo fastq
	fqname2=`basename ${fastq2} .fastq`
	
	fqbase=`basename ${fastq1} .fastq | sed 's/_R1//'`;

	echo " Processing ${fqbase} ..."

	echo "      Creating graphics pre analysis - ${fqname1} (FastQC) ..."
	fastqc -o ${fastqc_pre_out} -t ${threads} $fastq1 > ${fastqc_pre_out}/${fqname1}.fastqc.out.txt  2> ${fastqc_pre_out}/${fqname1}.fastqc.err.txt

	if [ -e ${fastq2} ]; then
		echo "      Creating graphics pre analysis - ${fqname2} (FastQC) ..."
		fastqc -o ${fastqc_pre_out} -t ${threads} $fastq2 > ${fastqc_pre_out}/${fqname2}.fastqc.out.txt  2> ${fastqc_pre_out}/${fqname2}.fastqc.err.txt
	fi

	####################################
	#Trimandos adaptadores e primers com scythe e cutadapt
	
	# scythe

	echo "      3p Adapters trimming for sample: ${fqname1} (Scythe) ..."
        scythe -p 0.3 -m ${scythe_out}/${fqname1}.scythe.matches.txt -o ${scythe_out}/${fqname1}.scythe.fastq -n 5 -M 1 -q sanger -a ${scythe_adapters_3p_R1_path} -i ${fastq1} 1> ${scythe_out}/${fqname1}.scythe.out.txt 2> ${scythe_out}/${fqname1}.scythe.err.txt
	if [ -e ${fastq2} ]; then
		echo "      3p Adapters trimming for sample: ${fqname2} (Scythe) ..."
	        scythe -p 0.3 -m ${scythe_out}/${fqname2}.scythe.matches.txt -o ${scythe_out}/${fqname2}.scythe.fastq -n 5 -M 1 -q sanger -a ${scythe_adapters_3p_R2_path} -i ${fastq2} 1> ${scythe_out}/${fqname2}.scythe.out.txt 2> ${scythe_out}/${fqname2}.scythe.err.txt
	fi

	# CutAdapt 
	
	cutadapt_adapters_5p_R1=`cat ${cutadapt_adapters_5p_R1_path} | perl -lane 'INIT {our @l=(); } if (($_)&&($_!~/^>/)) { push(@l, "--front=$_"); } END { print join(" ", @l);} '`
	cutadapt_adapters_5p_R2=`cat ${cutadapt_adapters_5p_R2_path} | perl -lane 'INIT {our @l=(); } if (($_)&&($_!~/^>/)) { push(@l, "--front=$_"); } END { print join(" ", @l);} '`
	cutadapt_fragments=`cat ${scythe_adapters_3p_R1_path} ${cutadapt_adapters_5p_R1_path} ${scythe_adapters_3p_R2_path} ${cutadapt_adapters_5p_R2_path} | perl -lane 'INIT {our @l=(); } if (($_)&&($_!~/^>/)) { push(@l, "--anywhere=$_"); } END { print join(" ", @l);} '`

	echo "      5p Adapter trimming for sample: ${fqname1} (CutAdapt) ..."
        cutadapt --format=fastq ${cutadapt_adapters_5p} --error-rate=0.1 --times=2 --overlap=4 --minimum-length=50 --output=${cutadapt_out}/${fqname1}.scythe.cutadapt5p.fastq ${scythe_out}/${fqname1}.scythe.fastq 1> ${cutadapt_out}/${fqname1}.scythe.cutadapt5p.out.txt 2> ${cutadapt_out}/${fqname1}.scythe.cutadapt5p.err.txt
	
	if [ -e ${fastq2} ]; then
		echo "      5p Adapter trimming for sample: ${fqname2} (CutAdapt) ..."
	        cutadapt --format=fastq ${cutadapt_adapters_5p} --error-rate=0.1 --times=2 --overlap=4 --minimum-length=50 --output=${cutadapt_out}/${fqname2}.scythe.cutadapt5p.fastq ${scythe_out}/${fqname2}.scythe.fastq 1> ${cutadapt_out}/${fqname2}.scythe.cutadapt5p.out.txt 2> ${cutadapt_out}/${fqname2}.scythe.cutadapt5p.err.txt
	fi

	echo "      Filtering irregular fragment reads for sample: ${fqname1} (CutAdapt) ..."
	cutadapt --format=fastq ${cutadapt_fragments} --discard-trimmed --overlap=15 --minimum-length=50 --output=${cutadapt_out}/${fqname1}.scythe.cutadapt5p.filtered.fastq ${cutadapt_out}/${fqname1}.scythe.cutadapt5p.fastq 1> ${cutadapt_out}/${fqname1}.scythe.cutadapt5p.filtered.out.txt 2> ${cutadapt_out}/${fqname1}.scythe.cutadapt5p.filtered.err.txt
	
	if [ -e ${fastq2} ]; then
		echo "      Filtering irregular fragment reads for sample: ${fqname2} (CutAdapt) ..."
		cutadapt --format=fastq ${cutadapt_fragments} --discard-trimmed --overlap=15 --minimum-length=50 --output=${cutadapt_out}/${fqname2}.scythe.cutadapt5p.filtered.fastq ${cutadapt_out}/${fqname2}.scythe.cutadapt5p.fastq 1> ${cutadapt_out}/${fqname2}.scythe.cutadapt5p.filtered.out.txt 2> ${cutadapt_out}/${fqname2}.scythe.cutadapt5p.filtered.err.txt
	fi

	echo "      Quality trimming for sample: ${fqbase} (PrinSeq) ..."
	
	fastqin="-fastq ${cutadapt_out}/${fqname1}.scythe.cutadapt5p.filtered.fastq"
	if [ -e ${fastq2} ]; then
		fastqin="${fastqin} -fastq2 ${cutadapt_out}/${fqname2}.scythe.cutadapt5p.filtered.fastq"
	fi

	prinseq-lite.pl -verbose ${fastqin} -out_format 3 -out_good ${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq \
			-out_bad null -qual_noscale -min_len 50 -ns_max_p 25 -noniupac -trim_tail_left 5 -trim_tail_right 5 -trim_qual_right 25 \
			-trim_qual_type mean -trim_qual_rule lt -trim_qual_window 3 -trim_qual_step 1 -lc_method dust -lc_threshold 30 \
		       	-no_qual_header 1 \
			>  ${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq.out.txt \
			2> ${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq.err.txt

	
	echo "      Creating graphics pos analysis - ${fqbase} (FastQC) ..."
	if [ ! -e ${fastq2} ]; then
		fastqc -o ${fastqc_pos_out} -t ${threads} ${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq.fastq > ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq.fastqc.out.txt  2> ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq.fastqc.err.txt
		
	else
		fastqc -o ${fastqc_pos_out} -t ${threads} ${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_1.fastq > ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_1.fastqc.out.txt  2> ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_1.fastqc.err.txt
		fastqc -o ${fastqc_pos_out} -t ${threads} ${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_2.fastq > ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_2.fastqc.out.txt  2> ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_2.fastqc.err.txt
	
		if [ -e "${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_1_singletons.fastq" ]; then
			fastqc -o ${fastqc_pos_out} -t ${threads} ${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_1_singletons.fastq > ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_1_singletons.fastqc.out.txt  2> ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_1_singletons.fastqc.err.txt
		fi
		if [ -e "${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_2_singletons.fastq" ]; then
			fastqc -o ${fastqc_pos_out} -t ${threads} ${prinseq_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_2_singletons.fastq > ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_2_singletons.fastqc.out.txt  2> ${fastqc_pos_out}/${fqbase}.scythe.cutadapt5p.filtered.prinseq_2_singletons.fastqc.err.txt
		fi
	fi

done

