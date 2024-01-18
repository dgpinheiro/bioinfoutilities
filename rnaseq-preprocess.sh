#!/bin/bash

indir=$1

# 1>&2 (redireciona STDOUT para STDERR)
# 2>&1 (redireciona STDERR para STDOUT)

# SE ${indir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO NA LINHA DE COMANDO
if [ ! ${indir} ]; then
	echo "[ERROR] Missing input directory." 1>&2
	exit
fi

# SE ${indir} NÃO É DIRETÓRIO
if [ ! -d ${indir} ]; then
	echo "[ERROR] Wrong input directory (${indir})." 1>&2
	exit
fi

outdir=$2

# SE ${outdir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO NA LINHA DE COMANDO
if [ ! ${outdir} ]; then
	echo "[ERROR] Missing output directory." 1>&2
	exit
fi

# SE ${outdir} NÃO É DIRETÓRIO
if [ ! -d ${outdir} ]; then
	echo "[ERROR] Wrong output directory (${outdir})." 1>&2
	exit
fi

# Número de CORES para o processamento
# ATENÇÃO: Não exceder o limite da máquina
THREADS=$3

if [ ! ${THREADS} ]; then
	echo "[ERROR] Missing number of threads." 1>&2
	exit
fi

if [[ ! ${THREADS} =~ [0-9]+ ]]; then
	echo "[ERROR] Wrong number of threads (${THREADS})." 1>&2
	exit
fi

if [ ! ${THREADS} -gt 0 ]; then
	echo "[ERROR] The number of threads must be greater than zero." 1>&2
	exit
fi

# Contaminantes
contaminants=$4

mkdir -p ${outdir}/processed/cleaned

if [ ${contaminants} ]; then

	if [ ! -e ${contaminants} ]; then
		echo "[ERROR] Wrong contaminants file (${contaminants})." 1>&2
		exit
	fi

	mkdir -p ${outdir}/processed/cleaned/refs
	
	curdir=`pwd`
	contaminants_abs_path=`readlink -f ${contaminants}`

	cd ${outdir}/processed/cleaned/refs

	if [ ! -e "./contaminants.fa" ]; then
		ln -s ${contaminants_abs_path} contaminants.fa
		if [ ! -e "./contaminants.1.bt2" ] ||
		   [ ! -s "./contaminants.1.bt2" ]; then
			echo -e "Creating bowtie2 index for contaminants (${contaminants}) ...\n"
			bowtie2-build	--threads ${THREADS} \
					contaminants.fa contaminants > bowtie2-build.log.out.txt 2> bowtie2-build.log.err.txt
		fi			
	fi

	cd ${curdir}
fi

# Criar estrutura de diretórios
mkdir -p ${outdir}/processed/fastqc/pre
mkdir -p ${outdir}/processed/atropos
mkdir -p ${outdir}/processed/prinseq
mkdir -p ${outdir}/processed/fastqc/pos

# É IMPORTANTE R1 E R2 TENHAM EXATAMENTE O MESMO NOME EXCETO ESTA PARTE: "1" ou "2" depois do "_R"
for r1 in `ls ${indir}/*_R1.fastq ${indir}/*_R1.fastq.gz 2>/dev/null`; do

	# R2 (caso sequenciamento paired-end)
	r2=`echo ${r1} | sed 's/_R1.fastq/_R2.fastq/'`

	# Este teste irá verificar a existência de R2
	if [ ! -e ${r2} ]; then
		echo "[ERROR] Not found mate for R1 (${r2})." 1>&2
		exit
	fi

	# Este teste irá verificar se os arquivos linkados no input não estão vazios
	if [ ! -s $(readlink -f ${r1}) ] || [ ! -s $(readlink -f ${r2}) ]; then
		echo "[ERROR] R1 or R2 input files were empty." 1>&2
		exit
	fi

	# OBTÉM O NOME (BASE) DA BIBLIOTECA EXCETO _R1.fastq ou R1.fastq.gz
        name=`basename ${r1} | sed 's/_R1.fastq\(\.gz\)\?//'`
	
	echo -e "Processing library ${name} ..."
	

	if [ ! -e "${outdir}/processed/fastqc/pre/${name}_R1_fastqc.html" ] || 
	   [ ! -s "${outdir}/processed/fastqc/pre/${name}_R1_fastqc.html" ]; then
        	
		echo -e "\tFastQC evaluation using sample ${name}: ${r1} ...\n"
	        
		fastqc -t ${THREADS} \
        	        ${r1} \
	                -o ${outdir}/processed/fastqc/pre/ \
        	        1> ${outdir}/processed/fastqc/pre/${name}_R1.log.out.txt \
                	2> ${outdir}/processed/fastqc/pre/${name}_R1.log.err.txt
	fi

	if [ ! -e "${outdir}/processed/fastqc/pre/${name}_R2_fastqc.html" ] ||	
	   [ ! -s "${outdir}/processed/fastqc/pre/${name}_R2_fastqc.html" ]; then
        	
		echo -e "\tFastQC evaluation using sample ${name}: ${r2} ...\n"

		fastqc -t ${THREADS} \
			${r2} \
			-o ${outdir}/processed/fastqc/pre/ \
			1> ${outdir}/processed/fastqc/pre/${name}_R2.log.out.txt \
			2> ${outdir}/processed/fastqc/pre/${name}_R2.log.err.txt
	fi

	if [ ! -e "${outdir}/processed/atropos/${name}_R1.atropos_final.fastq" ] ||
	   [ ! -s "${outdir}/processed/atropos/${name}_R1.atropos_final.fastq" ]; then

		echo -e "\tRunning atropos (insert) for adapter trimming using sample ${name}: ${r1} & ${r2} ...\n"

		# este trecho de ativação do ambiente "atropos" usando o pyenv pode ser modificado
		# para adaptar ao sistema instalado no servidor
		# no servidor hammer utilizamos o "pyenv"
		# necessário inicializar o "pyenv" dentro desta sessão
		
		#eval "$(pyenv init -)"
		# ativando ambiente Python para execução do atropos
		#pyenv activate atropos
		
		# A execução no modo insert só faz sentido usando sequenciamento paired-end
		atropos trim --aligner insert \
			     --error-rate 0.1 \
			     --times 2 \
			     --minimum-length 1 \
			     --op-order GAWCQ \
			     --match-read-wildcards \
			     --overlap 3 \
			     --quality-cutoff 10 \
			     --threads ${THREADS} \
			     --correct-mismatches conservative \
			     --pair-filter any \
			     -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTGAAAAA \
			     -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNNNNNNNNNNNNNGTGGTCGCCGTATCATTAAAAAA \
			     -o ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
			     -p ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
			     -pe1 ${r1} \
			     -pe2 ${r2} \
			     --untrimmed-output ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \
			     --untrimmed-paired-output ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq \
			       > ${outdir}/processed/atropos/${name}.atropos.log.out.txt \
			      2> ${outdir}/processed/atropos/${name}.atropos.log.err.txt
		
		echo -e "\tRunning atropos (adapter) for adapter trimming using sample ${name}: ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq & ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq ...\n"
		atropos trim    --aligner adapter \
			-e 0.1 \
			-n 2 \
			-m 1 \
			--match-read-wildcards \
			-O 3 \
			-q 10 \
			--pair-filter any \
			-pe1 ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \
			-pe2 ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq \
                        -a  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTGAAAAA \
                        -g  TTTTTTAATGATACGGCGACCACNNNNNNNNNNNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT \
                        -A  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNNNNNNNNNNNNNGTGGTCGCCGTATCATTAAAAAA \
                        -G  TTTTTCAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
			-T  ${THREADS} \
			-o  ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq  \
			-p  ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \
			 >  ${outdir}/processed/atropos/${name}.atropos_adapter.log.out.txt \
			2>  ${outdir}/processed/atropos/${name}.atropos_adapter.log.err.txt

		echo -e "\tMerging atropos adapter trimming results using sample ${name}: ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq and ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq + ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq and ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq ...\n"
		
		cat       ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
			  ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq \
			> ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq

		cat       ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
			  ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \
			> ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq

		echo -e "\tRemoving useless atropos results ...\n"

		rm -f ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
		      ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq \
		      ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \
		      ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
		      ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \
		      ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq
	fi

	if [ ! -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq" ] ||
	   [ ! -s "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq" ]	; then

		echo -e "\tPrinSeq processing: ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq & ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq ...\n"

		prinseq-lite.pl -fastq  ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq \
				-fastq2 ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq \
				-out_format 3 \
				-trim_qual_window 3 \
				-trim_qual_step 1 \
				-trim_qual_type mean \
				-trim_qual_rule lt \
				-trim_qual_right 30 \
				-out_good ${outdir}/processed/prinseq/${name}.atropos_final.prinseq \
				-out_bad  null \
				-lc_method dust \
				-lc_threshold 50 \
				-min_len 20 \
				-trim_tail_right 5 \
				-trim_tail_left 5 \
				-ns_max_p 80 \
				-min_qual_mean 28 \
				 > ${outdir}/processed/prinseq/${name}.atropos_final.prinseq.out.log \
				2> ${outdir}/processed/prinseq/${name}.atropos_final.prinseq.err.log
	
	fi

	if [ ! -e "${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_fastqc.zip" ] ||
	   [ ! -s "${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_fastqc.zip" ]	; then
		echo -e "\tFastQC pos-evaluation using sample ${name}: ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq\n"

		fastqc -t ${THREADS} \
		   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq \
		   -o ${outdir}/processed/fastqc/pos/ \
		    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1.log.out.txt \
		   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1.log.err.txt
	fi

	if [ ! -e "${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_fastqc.zip" ] ||
	   [ ! -s "${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_fastqc.zip" ]	; then
		echo -e "\tFastQC pos-evaluation using sample ${name}: ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2.fastq\n"

		fastqc -t ${THREADS} \
		   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2.fastq \
		   -o ${outdir}/processed/fastqc/pos/ \
		    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2.log.out.txt \
		   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2.log.err.txt
	fi
	
	# Remover (sempre) para garantir que não serão acrescentadas reads
	rm -f ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_singletons.fastq

	# SE EXISTIR <SAMPLE_NAME>.atropos_final.prinseq_1_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq" ]; then

		if [ ! -e "${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons_fastqc.zip" ] ||
		   [ ! -s "${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons_fastqc.zip" ]; then

			fastqc -t ${THREADS} \
			   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq \
			   -o ${outdir}/processed/fastqc/pos/ \
			    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons.log.out.txt \
			   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons.log.err.txt
		fi
		
		# concatena se existe "*_1_singletons.fastq"
		cat ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq > ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_singletons.fastq
	fi
	
	# SE EXISTIR <SAMPLE_NAME>.atropos_final.prinseq_2_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq" ]; then

		if [ ! -e "${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons_fastqc.zip" ] ||
		   [ ! -s "${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons_fastqc.zip" ]; then

			fastqc -t ${THREADS} \
			   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq \
			   -o ${outdir}/processed/fastqc/pos/ \
	        	    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_singletons.log.out.txt \
		           2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_singletons.log.err.txt
		fi
		
		# concatena se existe "*_2_singletons.fastq"
		cat ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq >> ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_singletons.fastq
	fi
		
	if [ ! -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq" ]; then
		touch ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq
	fi

	if [ ! -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq" ]; then
		touch ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq
	fi
	
	if [ ${contaminants} ]; then
	
		if [ ! -e "${outdir}/processed/cleaned/${name}.bowtie2.log.out.txt" ]; then

			bowtie2 -x ${outdir}/processed/cleaned/refs/contaminants \
				-1 ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq \
				-2 ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2.fastq \
				-U ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_singletons.fastq \
				1> ${outdir}/processed/cleaned/${name}.bowtie2.log.out.txt \
				2> ${outdir}/processed/cleaned/${name}.bowtie2.log.err.txt
			
			cat ${outdir}/processed/cleaned/${name}.bowtie2.log.out.txt | samtools view -F 4 - | cut -f 1 | sed 's/\/[12]$//' | sort | uniq | perl -lane 'print "$_/1"; print "$_/2"; print "$_";' > ${outdir}/processed/cleaned/${name}_exclude.txt

			rm -f 	${outdir}/processed/cleaned/${name}.atropos_final.prinseq.cleaned_1.fastq \
				${outdir}/processed/cleaned/${name}.atropos_final.prinseq.cleaned_2.fastq \
				${outdir}/processed/cleaned/${name}.atropos_final.prinseq.cleaned_1_singletons.fastq \
				${outdir}/processed/cleaned/${name}.atropos_final.prinseq.cleaned_2_singletons.fastq

			pullseq -i ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq \
				-n ${outdir}/processed/cleaned/${name}_exclude.txt -e \
				> ${outdir}/processed/cleaned/${name}.atropos_final.prinseq.cleaned_1.fastq
			
			pullseq -i ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2.fastq \
				-n ${outdir}/processed/cleaned/${name}_exclude.txt -e \
				> ${outdir}/processed/cleaned/${name}.atropos_final.prinseq.cleaned_2.fastq

			pullseq -i ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq \
				-n ${outdir}/processed/cleaned/${name}_exclude.txt -e \
				> ${outdir}/processed/cleaned/${name}.atropos_final.prinseq.cleaned_1_singletons.fastq
			
			pullseq -i ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq \
				-n ${outdir}/processed/cleaned/${name}_exclude.txt -e \
				> ${outdir}/processed/cleaned/${name}.atropos_final.prinseq.cleaned_2_singletons.fastq
		fi
	else
		curdir=`pwd`

		r1_cleaned=`readlink -f ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq`

		r2_cleaned=`readlink -f ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2.fastq`

		s1_cleaned=`readlink -f ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq`
		
		s2_cleaned=`readlink -f ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq`

		cd ${outdir}/processed/cleaned
			
		rm -f 	./${name}_exclude.txt \
			./${name}.bowtie2.log.out.txt \
			./${name}.bowtie2.log.err.txt \
			./${name}.atropos_final.prinseq.cleaned_1.fastq \
			./${name}.atropos_final.prinseq.cleaned_2.fastq \
			./${name}.atropos_final.prinseq.cleaned_1_singletons.fastq \
			./${name}.atropos_final.prinseq.cleaned_2_singletons.fastq
		
		ln -s ${r1_cleaned} ${name}.atropos_final.prinseq.cleaned_1.fastq
		ln -s ${r2_cleaned} ${name}.atropos_final.prinseq.cleaned_2.fastq
		ln -s ${s1_cleaned} ${name}.atropos_final.prinseq.cleaned_1_singletons.fastq
		ln -s ${s2_cleaned} ${name}.atropos_final.prinseq.cleaned_2_singletons.fastq

		cd ${curdir}
	fi

done



