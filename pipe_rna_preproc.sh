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

# as linhas que iniciam com cerquilha são comentários

#############################
# Caminhos dos arquivos e diretórios de entrada (input - diretório com os arquivos .fastq)
# 

# input - diretório contendo os arquivos de entrada no formato .fastq
input="$1"
if [ ! ${input} ]
then   
	echo "Missing data path to directory containing *.fastq files"
        exit
else   
        if [ ! -d ${input} ]
        then   
		echo "Wrong data path to input directory (${input}). Check first argument (Directory containing *.fastq files)"
                exit
        fi
fi


# Número de threads
threads=14

#############################
# Poda de adaptadores e primers
#

# Diretório base onde estão os adaptadores
adap_barc_primer_path="./refs/adapter/"

# scythe_adapters_3p_techseq_path - caminho para o arquivo contendo os adaptadores 3p
# (de acordo com o protocolo/tecnologia de sequenciamento) para o programa scythe
scythe_adapters_3p_techseq_path="${adap_barc_primer_path}/truseq_illumina_3p_adapt_R1.fa"

# cutadapt_adapters_5p_techseq_path - caminho para o arquivo contendo os adaptadores 3p 
# (de acordo com o protocolo/tecnologia de sequenciamento) para o programa scythe
cutadapt_adapters_5p_techseq_path="${adap_barc_primer_path}/truseq_illumina_5p_adapt_R1.fa"

#########################################
###
# Arquivos e diretórios de saída (output) 
#

# diretório a partir de onde o script será executado
base_out="."

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

echo "* Creating directories"

# Criando diretórios para as saídas dos programas que serão utilizados a seguir
mkdir -p ${fastqc_pre_out}
mkdir -p ${scythe_out}
mkdir -p ${cutadapt_out}
mkdir -p ${prinseq_out}
mkdir -p ${fastqc_pos_out}


###################################################
#FastQC
## percorrendo múltiplos valores retornados por um comando no shell (o comando irá listar todos os arquivos com extensão .fastq do diretório 1_raw

for fastq in `ls ${input}/*.fastq`; do


	# fqname - string com o nome do arquivo fastq
	fqname=`basename ${fastq} .fastq`

	echo " Processing ${fqname}..."

	echo "      Creating graphics pre analysis - ${fqname} (FastQC) ..."

	fastqc -o ${fastqc_pre_out} -t ${threads} $fastq > ${fastqc_pre_out}/${fqname}.fastqc.out.txt  2> ${fastqc_pre_out}/${fqname}.fastqc.err.txt

	####################################
	#Trimandos adaptadores e primers com scythe e cutadapt
		
	echo "      3p Adapters trimming for sample ${fqname} (scythe) ..."

        scythe -p 0.1 -m ${scythe_out}/${fqname}.scythe.matches.txt -o ${scythe_out}/${fqname}.scythe.fastq -n 5 -M 50 -q sanger -a ${scythe_adapters_3p_techseq_path} -i ${fastq} 1> ${scythe_out}/${fqname}.scythe.out.txt 2> ${scythe_out}/${fqname}.scythe.err.txt
	

	echo "      5p Adapter trimming for sample ${fqname} (cutadapt) ..."

	# cutadapt 
	cutadapt_adapters_5p=`cat ${cutadapt_adapters_5p_techseq_path} | perl -lane 'INIT {our @l=(); } if (($_)&&($_!~/^>/)) { push(@l, "--front=$_"); } END { print join(" ", @l);} '`
	
        cutadapt --format=fastq ${cutadapt_adapters_5p} --error-rate=0.1 --times=2 --overlap=4 --minimum-length=50 --output=${cutadapt_out}/${fqname}.scythe.cutadapt5p.fastq ${scythe_out}/${fqname}.scythe.fastq 1> ${cutadapt_out}/${fqname}.scythe.cutadapt5p.out.txt 2> ${cutadapt_out}/${fqname}.scythe.cutadapt5p.err.txt
	

	echo "      Filtering irregular fragment reads (cutadapt) ..."

	# cutadapt 
	cutadapt_fragments=`cat ${scythe_adapters_3p_techseq_path} ${cutadapt_adapters_5p_techseq_path} | perl -lane 'INIT {our @l=(); } if (($_)&&($_!~/^>/)) { push(@l, "--anywhere=$_"); } END { print join(" ", @l);} '`
	
	cutadapt --format=fastq ${cutadapt_fragments} --discard-trimmed --overlap=15 --minimum-length=50 --output=${cutadapt_out}/${fqname}.scythe.cutadapt5p.filtered.fastq ${cutadapt_out}/${fqname}.scythe.cutadapt5p.fastq 1> ${cutadapt_out}/${fqname}.scythe.cutadapt5p.filtered.out.txt 2> ${cutadapt_out}/${fqname}.scythe.cutadapt5p.filtered.err.txt

	echo "      Quality trimming for sample ${fqname} (prinseq) ..."
	
	prinseq-lite.pl -verbose -fastq ${cutadapt_out}/${fqname}.scythe.cutadapt5p.filtered.fastq -out_format 3 -out_good ${prinseq_out}/${fqname}.scythe.cutadapt5p.filtered.prinseq -out_bad null -qual_noscale -min_len 50 -ns_max_p 25 -noniupac -trim_tail_left 5 -trim_tail_right 5 -trim_qual_right 25 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 3 -trim_qual_step 1 -lc_method dust -lc_threshold 30  -no_qual_header 1> ${prinseq_out}/${fqname}.scythe.cutadapt5p.filtered.prinseq.out.txt 2> ${prinseq_out}/${fqname}.scythe.cutadapt5p.filtered.prinseq.err.txt
	
	echo "      Creating graphics pos analysis - ${fqname} (FastQC) ..."

	fastqc -o ${fastqc_pos_out} -t ${threads} ${prinseq_out}/${fqname}.scythe.cutadapt5p.filtered.prinseq.fastq > ${fastqc_pos_out}/${fqname}.fastqc.out.txt  2> ${fastqc_pos_out}/${fqname}.fastqc.err.txt

done;

