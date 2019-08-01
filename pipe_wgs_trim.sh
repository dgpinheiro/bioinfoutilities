#!/bin/bash

# Recebendo valores das variáveis a partir dos argumentos 
# da linha de comando
r1file=$1
r2file=$2
outdir=$3

# Número de threads (subprocessos) que serão lançados em paralelo
num_threads=2

# Testes para avaliar se parâmetros foram passados e 
# se arquivos (-e) ou diretórios (-d) existem
# Testes para avaliar se os nomes dos arquivos coincidem com
# os sufixos esperados "_R1.fastq" e "_R2.fastq"
if [ ! ${r1file} ]; then
	echo "Missing R1 file"
	exit
else
	if [ ! -e ${r1file} ]; then
		echo "Wrong R1 file (${r1file})"
		exit
	else
		if [[ ! ${r1file} =~ _R1.fastq$ ]]; then
			echo "Wrong name for R1 file (${r1file})"
			exit
		fi
	fi
fi


if [ ! ${r2file} ]; then
	echo "Missing R2 file"
	exit
else
	if [ ! -e ${r2file} ]; then
		echo "Wrong R2 file (${r2file})"
		exit
	else
		if [[ ! ${r2file} =~ _R2.fastq$ ]]; then
		echo "Wrong name for R2 file (${r2file})"
			exit
		fi
	fi
fi


if [ ! ${outdir} ]; then
	echo "Missing output directory"
	exit
else
	if [ ! -d ${outdir} ]; then
		echo "Wrong output directory"
		exit
	fi
fi

# Capturando os nomes dos arquivos, ou seja, sem os diretórios, e, também, 
# sem a extensão .fastq
bn1=`basename ${r1file} .fastq`
bn2=`basename ${r2file} .fastq`
# Substituindo no nome do arquivo R1 o trecho "_R1" por nada para ser
# o nome base independente de ser R1 ou R2, por exemplo, no caso dos arquivos
# de LOG para atropos e prinseq
bn=`echo ${bn1} | sed 's/_R1//'`

echo "R1: ${r1file}  [${bn1}]"
echo "R2: ${r2file}  [${bn2}]"
echo "OUTPUT DIRECTORY: ${outdir}"

# Contagem dos dados brutos
raw_countr1=`echo "$(cat ${r1file} | wc -l )/4" | bc`
raw_countr2=`echo "$(cat ${r2file} | wc -l )/4" | bc`

echo "R1 COUNT (raw): ${raw_countr1}"
echo "R2 COUNT (raw): ${raw_countr2}"
echo -e "${bn}_R1 COUNT (raw):${raw_countr1}\n${bn}_R2 COUNT (raw): ${raw_countr2}" > ${outdir}/${bn}-raw_count.txt

# Processamento da análise de qualidade com FastQC

mkdir -p ${outdir}/processed/fastqc/pre

fastqc ${r1file} -t ${num_threads} \
		 -o ${outdir}/processed/fastqc/pre/ \
		 &> ${outdir}/processed/fastqc/pre/${bn1}.fastqc.log

fastqc ${r2file} -t ${num_threads} \
		 -o ${outdir}/processed/fastqc/pre/ \
		 &> ${outdir}/processed/fastqc/pre/${bn2}.fastqc.log

# Processamento para retirada de sequências de adaptadores

mkdir -p ${outdir}/processed/atropos

# Inicialização do ambiente python apropriado para o atropos (hammer)
#eval "$(pyenv init -)"
#pyenv activate atropos

# source activate atropos

# Primeira análise com atropos utilizando método de alinhamento dos insertos 
# (sequência de interesse) de R1 e R2
atropos trim 	-pe1 ${r1file} \
		-pe2 ${r2file} \
		-a CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
		-A CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
		--aligner insert -e 0.1 -n 2 -m 20 --op-order GAWCQ\
	        --match-read-wildcards -O 20 -q 20 --pair-filter any \
		-T ${num_threads} --correct-mismatches conservative \
		-o ${outdir}/processed/atropos/${bn1}.atropos_insert.fastq  \
		-p ${outdir}/processed/atropos/${bn2}.atropos_insert.fastq  \
		--untrimmed-output ${outdir}/processed/atropos/${bn1}.atropos_untrimmed.fastq \
		--untrimmed-paired-output ${outdir}/processed/atropos/${bn2}.atropos_untrimmed.fastq \
		&> ${outdir}/processed/atropos/${bn}.atropos_insert.log.txt

# Compactando os arquivos com as leituras brutas
gzip -f ${r1file}
gzip -f ${r2file}

# Criando arquivos vazios (sequências que não foram trimadas - untrimmed) caso 
# ainda não tenham sido criados pelo atropos
if [ ! -e "${outdir}/processed/atropos/${bn1}.atropos_untrimmed.fastq" ]; then
	touch ${outdir}/processed/atropos/${bn1}.atropos_untrimmed.fastq
fi
if [ ! -e "${outdir}/processed/atropos/${bn2}.atropos_untrimmed.fastq" ]; then
	touch ${outdir}/processed/atropos/${bn2}.atropos_untrimmed.fastq
fi


# Segunda análise com atropos utilizando método de alinhamento com adaptadores,
# somente para os arquivos não trimados
atropos trim 	-pe1 ${outdir}/processed/atropos/${bn1}.atropos_untrimmed.fastq \
		-pe2 ${outdir}/processed/atropos/${bn2}.atropos_untrimmed.fastq \
		-a CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
		-g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
		-A CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
		-T ${num_threads} \
		--aligner adapter -e 0.1 -n 2 -m 20 \
		--match-read-wildcards -O 3 -q 20 --pair-filter both \
		-o ${outdir}/processed/atropos/${bn1}.atropos_adapter.fastq  \
		-p ${outdir}/processed/atropos/${bn2}.atropos_adapter.fastq \
		&> ${outdir}/processed/atropos/${bn}.atropos_adapter.log.txt

# Concatenando arquivos resultantes da primeira e segunda análise com atropos
# para um único arquivo "*.atropos_final.fastq"
cat 	${outdir}/processed/atropos/${bn1}.atropos_insert.fastq \
	${outdir}/processed/atropos/${bn1}.atropos_adapter.fastq \
	> ${outdir}/processed/atropos/${bn1}.atropos_final.fastq

# Removendo arquivos já concatenados em *_R1.atropos_final.fastq e o _R1.atropos_untrimmed.fastq
#rm -f 	${outdir}/processed/atropos/${bn1}.atropos_insert.fastq \
#	${outdir}/processed/atropos/${bn1}.atropos_adapter.fastq \
#	${outdir}/processed/atropos/${bn1}.atropos_untrimmed.fastq

cat	${outdir}/processed/atropos/${bn2}.atropos_insert.fastq \
	${outdir}/processed/atropos/${bn2}.atropos_adapter.fastq \
	> ${outdir}/processed/atropos/${bn2}.atropos_final.fastq

# Removendo arquivos já concatenados em *_R2.atropos_final.fastq
#rm -f	${outdir}/processed/atropos/${bn2}.atropos_insert.fastq \
#	${outdir}/processed/atropos/${bn2}.atropos_adapter.fastq \
#	${outdir}/processed/atropos/${bn2}.atropos_untrimmed.fastq

# desativando ambiente python do atropos
#pyenv deactivate

# deactivate atropos

# Contagem de sequências pós análises com atropos
atropos_countr1=`echo "$(cat ${outdir}/processed/atropos/${bn1}.atropos_final.fastq | wc -l )/4" | bc`
atropos_countr2=`echo "$(cat ${outdir}/processed/atropos/${bn2}.atropos_final.fastq | wc -l )/4" | bc`

echo "R1 COUNT (atropos): ${atropos_countr1}"
echo "R2 COUNT (atropos): ${atropos_countr2}"
echo -e "${bn}_R1 COUNT (atropos):${atropos_countr1}\n${bn}_R2 COUNT (atropos): $atropos_countr2}" > ${outdir}/${bn}-atropos_count.txt

# Processamento para poda (trimagem) de regiões de baixa qualidade, cauda Poli-A/T e filtragem
# de sequências pequenas e de baixa complexidade
prinseq_out="${outdir}/processed/prinseq"

mkdir -p ${prinseq_out}

prinseq-lite.pl         -verbose \
			-fastq ${outdir}/processed/atropos/${bn1}.atropos_final.fastq \
                        -fastq2 ${outdir}/processed/atropos/${bn2}.atropos_final.fastq \
			-min_len 50 \
			-out_format 3 \
			-trim_qual_step 1 \
			-trim_qual_window 10\
			-trim_qual_type min \
			-trim_qual_right 20 \
			-trim_qual_rule lt \
			-out_good ${prinseq_out}/${bn}.prinseq \
			-out_bad null \
			-lc_method dust \
			-lc_threshold 30 \
			1> ${prinseq_out}/${bn}.prinseq.log.out.txt \
			2> ${prinseq_out}/${bn}.prinseq.log.err.txt

# Compactar os arquivos já utilizados pelo prinseq-lite.pl
gzip -f ${outdir}/processed/atropos/${bn1}.atropos_final.fastq
gzip -f ${outdir}/processed/atropos/${bn2}.atropos_final.fastq

# Renomeando as saídas do prinseq para obedecerem o padrão estabelecido por
# convenção (disciplina)
mv ${outdir}/processed/prinseq/${bn}.prinseq_1.fastq ${outdir}/processed/prinseq/${bn1}.atropos.prinseq.fastq
mv ${outdir}/processed/prinseq/${bn}.prinseq_2.fastq ${outdir}/processed/prinseq/${bn2}.atropos.prinseq.fastq

# Criando os arquivos singletons caso não tenham sido gerados, o que serve 
# para evitar erros futuros provocados pela falta de algum desses arquivos
if [ ! -e ${outdir}/processed/prinseq/${bn}.prinseq_1_singletons.fastq ]; then
	touch ${outdir}/processed/prinseq/${bn}.prinseq_1_singletons.fastq
fi

if [ ! -e ${outdir}/processed/prinseq/${bn}.prinseq_2_singletons.fastq ]; then
	touch ${outdir}/processed/prinseq/${bn}.prinseq_2_singletons.fastq
fi

# Renomeando singletons
mv ${outdir}/processed/prinseq/${bn}.prinseq_1_singletons.fastq ${outdir}/processed/prinseq/${bn1}_singletons.atropos.prinseq.fastq
mv ${outdir}/processed/prinseq/${bn}.prinseq_2_singletons.fastq ${outdir}/processed/prinseq/${bn2}_singletons.atropos.prinseq.fastq

# Contagem das sequências após processamento com prinseq para paired-end (pe) e single-end (se)
prinseq_pe_countr1=`echo "$(cat ${outdir}/processed/prinseq/${bn1}.atropos.prinseq.fastq | wc -l )/4" | bc`
prinseq_pe_countr2=`echo "$(cat ${outdir}/processed/prinseq/${bn2}.atropos.prinseq.fastq | wc -l )/4" | bc`

prinseq_se_countr1=`echo "$(cat ${outdir}/processed/prinseq/${bn1}_singletons.atropos.prinseq.fastq | wc -l )/4" | bc`
prinseq_se_countr2=`echo "$(cat ${outdir}/processed/prinseq/${bn2}_singletons.atropos.prinseq.fastq | wc -l )/4" | bc`

echo "R1 COUNT PE (prinseq): ${prinseq_pe_countr1}"
echo "R2 COUNT PE (prinseq): ${prinseq_pe_countr2}"
echo -e "${bn}_R1 COUNT PE (prinseq): ${prinseq_pe_countr1}\n${bn}_R2 COUNT PE (prinseq): ${prinseq_pe_countr2}" > ${outdir}/${bn}-prinseq_pe_count.txt

echo "R1 COUNT SE (prinseq): ${prinseq_se_countr1}"
echo "R2 COUNT SE (prinseq): ${prinseq_se_countr2}"
echo -e "${bn}_R1 COUNT SE (prinseq): ${prinseq_se_countr1}\n${bn}_R2 COUNT SE (prinseq): ${prinseq_se_countr2}" > ${outdir}/${bn}-prinseq_se_count.txt

# Processamento com FastQC para análise de qualidade após etapas anteriores
# de pré-processamento (atropos e prinseq)
mkdir -p ${outdir}/processed/fastqc/pos

# Obtendo número aleatório para criação de diretório temporário (até análise com fastqc)
rand=${RANDOM}
mkdir -p /dev/shm/${rand}/

# Concatenando arquivos com sequências pareadas e não pareads (singletons)
cat ${outdir}/processed/prinseq/${bn1}.atropos.prinseq.fastq  ${outdir}/processed/prinseq/${bn1}_singletons.atropos.prinseq.fastq > /dev/shm/${rand}/${bn1}_all.atropos.prinseq.fastq
cat ${outdir}/processed/prinseq/${bn2}.atropos.prinseq.fastq  ${outdir}/processed/prinseq/${bn2}_singletons.atropos.prinseq.fastq > /dev/shm/${rand}/${bn2}_all.atropos.prinseq.fastq

fastqc 	/dev/shm/${rand}/${bn1}_all.atropos.prinseq.fastq \
	-t ${num_threads} \
	-o ${outdir}/processed/fastqc/pos/ \
	&> ${outdir}/processed/fastqc/pos/${bn1}.fastqc.log

fastqc 	/dev/shm/${rand}/${bn2}_all.atropos.prinseq.fastq \
	-t ${num_threads} \
	-o ${outdir}/processed/fastqc/pos/ \
	&> ${outdir}/processed/fastqc/pos/${bn2}.fastqc.log

# Removendo diretório temporário após análise com FastQC
rm -fr /dev/shm/${rand}/


# Final do arquivo

