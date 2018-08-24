#!/bin/bash

declare -A MYDESIGN

####
### Loading config file
####

if [ -e "$1" ]; then
	configfile=$1
	# execopt - align: somente alinhamento; diff: análise de expressão diferencial; all: todas
	execopt=$2
	# input - diretório contendo os arquivos de entrada no formato .fastq
	input=$3
	# aligner - alinhador [Default: tophat2]
	aligner=$4
	# basedir_out - diretório de saída
	basedir_out=$5
else
	configfile='./cuffpipe.cfg'
	# execopt - align: somente alinhamento; diff: análise de expressão diferencial; all: todas
	execopt=$1
	# input - diretório contendo os arquivos de entrada no formato .fastq
	input=$2
	# aligner - alinhador [Default: tophat2]
	aligner=$3
	# basedir_out - diretório de saída
	basedir_out=$4
fi

if [ ! -e ${configfile} ]; then
	echo "Not found configuration file (${configfile})"
	exit
else
	dos2unix ${configfile}
fi

echo "Using ${configfile}"

configfile_secured=$(mktemp)

# checagem se o arquivo contem algo indesejavel (que não segue o padrao esperado para arquivos de configuracao chave=valor)
if egrep -q -v '^[^ ]*=[^;]*' "${configfile}"; then
  echo "Config file is unclean, cleaning it..." >&2
  # filtragem do arquivo original para um novo arquivo
  egrep '^[^ ]*=[^;&]*'  "${configfile}" > "${configfile_secured}"
  configfile="${configfile_secured}"
fi

# carrega o arquivo original ou sua versao filtrada
source "${configfile}"

# remove arquivo de configuração temporario
rm ${configfile_secured}


if [ ! ${execopt} ]
then
	echo "Missing execution option (\"align\", \"diff\" or \"all\")"
	exit 1;
else
	if [ ! ${execopt} == "align" ] && [ ! ${execopt} == "diff" ] && [ ! ${execopt} == "all" ]
	then
		echo "Wrong option: ${execopt}"
		exit 1;
	fi
fi


# as linhas que iniciam com cerquilha são comentários

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

if [ ! ${aligner} ]
then
	aligner="tophat2"
	read -p "Press [Enter] key to start with ${aligner} or Ctrl+C to include other aligner in command line ..."
else
	if [ ! ${aligner} == "tophat2" ] && [ ! ${aligner} == "star" ] && [ ! ${aligner} == "bbmap" ] && [ ! ${aligner} == "hisat2" ]
	then
		echo "Wrong aligner (${aligner})"
		exit
	fi
fi




###
# Arquivos e diretórios de entrada (input)
#


# num_threads - inteiro com o número de processadores para serem utilizados nos parâmetros dos programas que serão executados a seguir
# Checando se num_threads foi definida
if [ ! ${num_threads} ]; then
	num_threads=1
	echo "Using default number of threads (${num_threads})";
fi

# refgff - caminho para o arquivo contendo as coordenadas gênicas
# Checando se refgff foi definida
if [ ! ${refgff} ]; then
	echo "Reference gff file (refgff) was not defined";
	exit;
fi
	
# check if cufflinks_library_type was defined
if [ ! ${cufflinks_library_type} ]; then
	cufflinks_library_type="fr-unstranded"
	echo "Cufflinks library type (cufflinks_library_type) was not defined. Using default: ${cufflinks_library_type}"
fi

if [ ${aligner} == "tophat2" ]; then
	# refbt2idx - caminho para a base com o índice do genoma para o bowtie2
	# Checando se refbt2idx foi definida
	if [ ! ${refbt2idx} ]; then
		echo "Reference Bowtie2 index (refbt2idx) was not defined";
		exit;
	fi

	# tophat_transcriptome_index - bowtie2 transcriptome index relacionado ao genes.gtf selecionado
	# Checando se tophat_transcriptome_index foi definida
	if [ ! ${tophat_transcriptome_index} ]; then
		tophat_transcriptome_index="./refs/transcriptome_index/transcriptome"
		echo "Using default tophat_transcriptome_index (${tophat_transcriptome_index})";
	fi
	mkdir -p $(dirname ${tophat_transcriptome_index})
	
else
	if [ ${aligner} == "star" ]; then
		# refstaridx - caminho para a base com o índice do genoma para o STAR
		# Checando se refstaridx foi definida
		if [ ! ${refstaridx} ]; then
			echo "Reference STAR index (refstaridx) was not defined";
			exit;
		fi
	else
		if [ ${aligner} == "bbmap" ]; then
			# refbbmapidx - caminho para a base com o índice do genoma para o bbmap
			# Checando se refbbmapidx foi definida
			if [ ! ${refbbmapidx} ]; then
				echo "Reference bbmap index (refbbmapidx) was not defined";
				exit;
			fi
		else
			if [ ${aligner} == "hisat2" ]; then
				if [ ! ${refhs2idx} ]; then
					echo "Reference hisat index (refhs2idx) was not defined";
					exit;
				fi
				
				if [ ! ${hisat_rna_strandness} ]; then
					hisat_rna_strandness=""
					echo "Hisat2 strandness parameter (hisat_rna_strandness) was not defined. Using default: ${hisat_rna_strandness}"
				fi

				
			fi
		fi
	fi
fi

# refgenomefa - caminho para o arquivo com as sequências do genoma no formato multifasta
# Checando se refgenomefa foi definida
if [ ! ${refgenomefa} ]; then
	echo "Reference genome fasta file (refgenomefa) was not defined";
	exit;
fi


# scythe_adapters_r1_path - caminho dos adaptadore R1 Scythe
# Checando se scythe_adapters_r1_path foi definida
if [ ! ${scythe_adapters_r1_path} ]; then
	echo "scythe_adapters_r1_path not defined";
	exit;
fi

# scythe_adapters_r2_path - caminho dos adaptadore R2 Scythe
# Checando se scythe_adapters_r2_path foi definida
if [ ! ${scythe_adapters_r2_path} ]; then
	echo "scythe_adapters_r2_path not defined";
	exit;
fi


# cutadapt_5p_adapters_r1_path - caminho dos adaptadore R1 CutAdapt
# Checando se cutadapt_5p_adapters_r1_path foi definida
if [ ! "${cutadapt_5p_adapters_r1_path}" ]; then
	echo "cutadapt_5p_adapters_r1_path not defined";
	exit;
fi

# cutadapt_5p_adapters_r2_path - caminho dos adaptadore R1 CutAdapt
# Checando se cutadapt_5p_adapters_r2_path foi definida
if [ ! "${cutadapt_5p_adapters_r1_path}" ]; then
	echo "cutadapt_5p_adapters_r2_path not defined";
	exit;
fi

if [ ! "${input_ext}" ]; then
	echo "input_ext not defined. Using DEFAULT: .fastq";
	input_ext=".fastq"
fi

if [ ! "${read_cmd}" ]; then
	echo "read_cmd not defined. Using DEFAULT: cat ";
fi


###
# Arquivos e diretórios de saída (output) 
#


if [ ! ${basedir_out} ]; then
	basedir_out="."
	read -p "Press [Enter] key to start with output path directory ${basedir_out} or Ctrl+C to include other path to output directory in command line ..."
else
	if [ ! -e ${basedir_out} ]; then
		mkdir -p ${$basedir_out}	
	fi		
fi

echo " * Using Output Directory ${basedir_out}"



# fastx_out - caminho para o diretório de saída do fastx_trimmer
fastx_out="${basedir_out}/raw2"

# cutadapt_out - caminho para o diretório de saída do cutadapt
cutadapt_out="${basedir_out}/processed/cutadapt"

# scythe_out - caminho para o diretório de saída do scythe
scythe_out="${basedir_out}/processed/scythe"

# seqpurge_out - caminho para o diretório de saída do SeqPurge
seqpurge_out="${basedir_out}/processed/seqpurge"

# align0_out - caminho para o diretório de saída do alinhamento (alignment 0)
align0_out="${basedir_out}/align0"
	
# align1_out - caminho para o diretório de saída do alinhamento (alignment 1)
align1_out="${basedir_out}/align1"

# align2_out - caminho para o diretório de saída do alinhamento (alignment 2)
align2_out="${basedir_out}/align2"

# align3_out - caminho para o diretório de saída do alinhamento (alignment 3)
align3_out="${basedir_out}/align3"

# align4_out - caminho para o diretório de saída do alinhamento (alignment 4)
align4_out="${basedir_out}/align4"

# finalalign_out - caminho para o diretório de saída do alinhamento final
finalalign_out="${basedir_out}/finalalign"

# prinseq_out - caminho para o diretório de saída do prinseq
prinseq_out="${basedir_out}/processed/prinseq"

refgff_format=`echo ${refgff} | perl -ne 'chomp; ~/\.([^\.]+)$/; print $1;'`

# biosamps - array com os nomes dos arquivos .fastq de cada réplica (técnica/biológica)
biosamps=()
# biogroups - array com os nomes dos grupos de amostras obtidos a partir do nome do arquivo fastq inicial
biogroups=()

# menor tamanho de sequência
min=1000

echo "* Identifying sequence size: "

lenvar=()
# percorrendo múltiplos valores retornados por um comando no shell (o comando irá listar todos os arquivos com extensão .fastq do diretório ${input})
for fastq in `ls ${input}/*R?${input_ext}`; do

	# obtendo nome do arquivo 
	fastqbn=`basename ${fastq} ${input_ext}`;

	# carregando valores de biosamps com o nome de cada réplica biológica. Por ex.: FLN_B1 obtido de FLN_B1_T1
	# Exemplo na linha de comando para remover valores duplicados * pois é possível que haja mais de 1 réplica técnica: 
	# LIST=("A1" "A1" "A1" "A2" "C2" "C1" "B1" "B2") && echo $(printf "%s\n" ${LIST[@]} | sort -u)
	biosamps=($( printf "%s\n" ${biosamps[@]} `echo ${fastqbn}  | cut -d '_' -f 1,2` | sort -u ))
	
	# carregando valores de biogroups com o nome de cada grupo de amostras. Por ex.: FLN obtido de FLN_B1_T1
	biogroups=($(printf "%s\n" ${biogroups[@]} `echo ${fastqbn}  | cut -d '_' -f 1` | sort -u ))
	
	# obtendo a primeira sequência do arquivo ${fastq}
	firstSeq=`${read_cmd} ${fastq} | head -2 |  tail -1`

	echo "   ${fastqbn}: ${#firstSeq} bases"

	# testando se o tamanho de ${firstSeq} é menor que mínimo
	if [ ${#firstSeq} -le ${min} ]
	then
		min=${#firstSeq}
	fi
	
	# carregando valores de tamanho de sequências
	lenvar=($(printf "%s\n" ${lenvar[@]} ${#firstSeq} | sort -n -u ))
done

echo "   Minimum sequence size: ${min}"

if [ ${execopt} == "align" ] || [ ${execopt} == "all" ]
then
	
	if [ ! ${max_intron_size} ] || [ ! ${min_intron_size} ]; then
		mkdir -p ${basedir_out}/introntab

		echo "* Gene structure statistics from current genome reference annotation (introntab.pl)"
		if [ ! -e "${basedir_out}/introntab/genome_annotation_stats.txt" ]; then
			echo "   Running introntab.pl based on ${refgff_format}"
			
			if [ ${refgff_format} == "gtf" ]; then
				gtfbn=`basename ${refgff} .gtf`
				if [ ! -e "${basedir_out}/introntab/${gtfbn}.gff" ]; then
					echo "      Converting gtf to gff for introntab.pl"
					gtf2gff3 ${refgff} 2> /dev/null > ${basedir_out}/introntab/${gtfbn}.gff 
				fi
				introntab.pl --format gff ${basedir_out}/introntab/${gtfbn}.gff > ${basedir_out}/introntab/genome_annotation_stats.txt
			else
				introntab.pl --format ${refgff_format} < ${refgff} > ${basedir_out}/introntab/genome_annotation_stats.txt
			fi
			
		fi

		max_intron_size=`cut -f 10 ${basedir_out}/introntab/genome_annotation_stats.txt | perl -lane 'next if (($.<=3)||($_=~/^#/)); my @isize=split(/,/, $_);  foreach my $s ( @isize ) { print $s; } ' | nsort -n | awk '{all[NR] = $0} END{print all[int(NR*0.75 - 0.5)]}' `
		min_intron_size=`cut -f 10 ${basedir_out}/introntab/genome_annotation_stats.txt | perl -lane 'next if (($.<=3)||($_=~/^#/)); my @isize=split(/,/, $_);  foreach my $s ( @isize ) { print $s; } ' | nsort -n | awk '{all[NR] = $0} END{print all[int(NR*0.25 - 0.5)]}' `

		if [ ! ${max_intron_size} ] || [ ! ${min_intron_size} ] ; then
			max_intron_size=`grep 'Genome_size:' ${basedir_out}/introntab/genome_annotation_stats.txt | perl -lane '$_=~/(\d+)/; print $1;' `;
			min_intron_size=`grep 'Genome_size:' ${basedir_out}/introntab/genome_annotation_stats.txt | perl -lane '$_=~/(\d+)/; print $1;' `;
			
			if [ ! ${max_intron_size} ] || [ ${min_intron_size} ]; then
				echo "   Not found intron size (PLEASE! Check introntab.pl output)"
				exit
			fi
		fi
		
	else
	
		if [ ${max_intron_size} ] || [ ${min_intron_size} ]; then
			echo "   Max intron size: ${max_intron_size}"
			echo "   Min intron size: ${min_intron_size}"
		else
			echo "   Not found intron size (PLEASE! Check min_intron_size and max_intron_size parameters)"
			exit
		fi
	fi
	
	echo "* Creating directories (alignments' step)"

	# Criando diretórios para as saídas dos programas que serão utilizados a seguir
	mkdir -p ${fastx_out}
	mkdir -p ${align0_out}
	mkdir -p ${align1_out}
	mkdir -p ${align2_out}
	mkdir -p ${align3_out}
	mkdir -p ${align4_out}
	mkdir -p ${finalalign_out}
	mkdir -p ${scythe_out}
	mkdir -p ${cutadapt_out}
	mkdir -p ${prinseq_out}
	
	
	for bs in ${biosamps[@]}; do
		# Concatena os dados das réplicas técnicas em uma única réplica biológica e podas as sequências para deixá-las do mesmo tamanho
		# caso ainda não tenham sido concatenados e podados (não existe -e)
		if [ ! -e "${fastx_out}/${bs}_R1.fastq.gz" ]; then
			if [ ! -e "${fastx_out}/${bs}_R1.fastq" ]  
			then
				techreps=(`ls ${input}/${bs}_T*_R1${input_ext}`)
				if ( ((${#techreps[@]} > 1)) || ((${#lenvar[@]} > 1)) )
				then
					if ((${#techreps[@]} > 1))
					then
						echo "* Trimming ${bs} by minimum length (${min}) and merging ${#techreps[@]} technical replicates of ${bs} [R1] ..."
					else
						echo "* Trimming ${bs} by minimum length (${min}) [R1] ..."
					fi
					rm -f ${fastx_out}/${bs}_R1.fastq
					for trinr1 in ${input}/${bs}_T*_R1${input_ext}; do
						trbn=`basename ${trinr1} _R1${input_ext}`
						${read_cmd} ${trinr1} | fastx_trimmer -l ${min} -Q 33 | awk -v TRBN="${trbn}" '{print (NR%4 == 1) ? "@" TRBN "_" ++i : ((NR%4 == 3) ? "+" : $0) }' >> ${fastx_out}/${bs}_R1.fastq
					done
				else
					link1=$(readlink -f ${input}/${bs}_T1_R1${input_ext})
					ftype=`file -b -i $link1 | cut -d ';' -f 1`

					if [ ${ftype} == "text/plain" ]; then
						echo "* Linking ${bs} [R1] ..."
						ln -s ${link1} ${fastx_out}/${bs}_R1.fastq
					else
						if [ ! -e "${fastx_out}/${bs}_R1.fastq.gz" ]; then
							echo "* Reading and writing ($ftype) ${bs} [R1] ..."
							${read_cmd} ${link1} > ${fastx_out}/${bs}_R1.fastq
						else
							echo "* Compressed ${bs} [R1]"
						fi
					fi
				fi
			fi
			if [ ! -e "${fastx_out}/${bs}_R2.fastq" ]  
			then
				# Erro esta sendo redirecionado para /dev/null pois pode não haver arquivos com o padrão desejado, o que nesse caso seriam 
				# bibliotecas single-end, ou seja, somente possuem "*_T*_R1.fastq"
				techreps=(`ls ${input}/${bs}_T*_R2${input_ext} 2>/dev/null`)
				if ((${#techreps[@]} > 0))
				then
					if ( ((${#techreps[@]} > 1)) || ((${#lenvar[@]} > 1)) )
					then
						if ((${#techreps[@]} > 1))
						then
							echo "* Trimming ${bs} by minimum length (${min}) and merging ${#techreps[@]} technical replicates of ${bs} [R2] ..."
						else
							echo "* Trimming ${bs} by minimum length (${min}) [R2] ..."
						fi
						for trinr2 in ${input}/${bs}_T*_R2${input_ext}; do
							trbn=`basename ${trinr2} _R2${input_ext}`
							${read_cmd} ${trinr2} | fastx_trimmer -l ${min} -Q 33 | awk -v TRBN="${trbn}" '{print (NR%4 == 1) ? "@" TRBN "_" ++i : ((NR%4 == 3) ? "+" : $0) }' >> ${fastx_out}/${bs}_R2.fastq
						done
					else
						link2=$(readlink -f ${input}/${bs}_T1_R2${input_ext})
						ftype=`file -b -i ${link2} | cut -d ';' -f 1`

						if [ ${ftype} == "text/plain" ]; then
							echo "* Linking ${bs} [R2] ..."
							ln -s ${link2} ${fastx_out}/${bs}_R2.fastq
						else
							if [ ! -e "${fastx_out}/${bs}_R2.fastq.gz" ]; then
								echo "* Reading and writing ($ftype) ${bs} [R2] ..."
								${read_cmd} ${link2} > ${fastx_out}/${bs}_R2.fastq
							else
								echo "* Compressed ${bs} [R2]"
							fi					
						fi
					fi
				fi
			fi
		fi
	done
	
	function hisat2_alignment_step {
		local samplename=$1
		local step=$2
		local alignX_sample_out=$3
		local alignL_sample_out=$4
		local fastqin=$5
		local fastq2in=$6
		local fastquin=$7
		
		#local nonovel=" --no-temp-splicesite "
		local nonovel=""
		local nomixed=" --no-mixed "
		local nodiscordant=" --no-discordant "
		local nosoftclip=" --no-softclip "
		if [ ${step} == 3 ]
		then
			nomixed=""
			nonovel=""
			nodiscordant=""
			nosoftclip=""
		fi
		
		mkdir -p ${alignX_sample_out}
				
		local rnastrandness=""
		if [ ${hisat_rna_strandness} ]; then
			rnastrandness="--rna-strandness ${hisat_rna_strandness}"
		fi
		if [ ! ${hisat_extra_params} ]; then
			hisat_extra_params=""
		fi

		local cmdhisat2="hisat2 ${rnastrandness} ${hisat_extra_params} --max-seeds 100 ${nomixed} ${nonovel} ${nodiscordant} ${nosoftclip} --threads=${num_threads} --dta-cufflinks --maxins 600 --min-intronlen  ${min_intron_size} --max-intronlen ${max_intron_size} --novel-splicesite-outfile ${alignX_sample_out}/junctions.junc "
		
		if [ -s "${alignL_sample_out}/all_junctions.junc" ]
		then   
			cmdhisat2="${cmdhisat2} --novel-splicesite-infile ${alignL_sample_out}/all_junctions.junc "
		fi
		
		cmdhisat2="${cmdhisat2} -x ${refhs2idx}"
		
		if [ -e ${fastq2in} ] 
		then
			cmdhisat2="${cmdhisat2} -1 ${fastqin} -2 ${fastq2in}"
			if [ ${fastquin} ]; then
				cmdhisat2="${cmdhisat2} -U ${fastquin}"
			fi
		else
			cmdhisat2="${cmdhisat2} -U ${fastqin}"
		fi
		
		cmdhisat2="${cmdhisat2} 2>${alignX_sample_out}/${samplename}.hisat2.log | samtools view -@ ${num_threads} -F 4 -b - > ${alignX_sample_out}/accepted_hits.bam 2> ${alignX_sample_out}/${samplename}.samtools.log.txt"
		
		#echo ${cmdhisat2}
		eval ${cmdhisat2}
		
		validaln=$(samtools view ${alignX_sample_out}/accepted_hits.bam  | head -1)
		
		if [ ! "${validaln}" ]; then
			echo "Warning: not found any valid alignment (${alignX_sample_out}/accepted_hits.bam) !"
		fi
		
		if [ -e "${alignX_sample_out}/junctions.junc" ]; then
			cat ${alignX_sample_out}/junctions.junc | nsort -u > ${alignX_sample_out}/all_junctions.junc
		fi
		
		samtools view -@ ${num_threads} -F 4 ${alignX_sample_out}/accepted_hits.bam | cut -f 1 | nsort -u > ${alignX_sample_out}/${samplename}.mapped.txt
		
		cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fastqin} > ${alignX_sample_out}/${samplename}.unmapped_R1.fastq
		
		if [ -e ${fastq2in} ]
		then
			cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fastq2in} > ${alignX_sample_out}/${samplename}.unmapped_R2.fastq
			if [ ${fastquin} ]
			then
				rm -f ${alignX_sample_out}/${samplename}.unmapped_S*.fastq
				rm -f ${alignX_sample_out}/${samplename}.unmapped_U.fastq
				for fquin in $(echo "${fastquin}" | tr "," "\n"); do
					if [ -e ${fquin} ]; then
						if [[ ${fquin} =~ "prinseq_1_singletons" ]] || [[ ${fquin} =~ "singleton_R1" ]]; then
							cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fquin} >> ${alignX_sample_out}/${samplename}.unmapped_S1.fastq
						elif [[ ${fquin} =~ "prinseq_2_singletons" ]] || [[ ${fquin} =~ "singleton_R2" ]]; then
							cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fquin} >> ${alignX_sample_out}/${samplename}.unmapped_S2.fastq
						else
							cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fquin} >> ${alignX_sample_out}/${samplename}.unmapped_U.fastq
						fi
					fi
				done
			fi
		fi
		
	}
	
	function tophat_alignment_step {
		local samplename=$1
		local step=$2
		local alignX_sample_out=$3
		local alignL_sample_out=$4
		local fastqin=$5
		local fastq2in=$6
		
		local nonovel=" --no-novel-juncs "
		local nomixed=" --no-mixed "
		local nodiscordant=" --no-discordant "
		if [ ${step} == 3 ]
		then
			nomixed=""
			nonovel=""
			nodiscordant=""
		fi
		
		local cmdtophat2="tophat2 ${nomixed} ${nonovel} ${nodiscordant} --num-threads=${num_threads} --library-type=${cufflinks_library_type} --prefilter-multihits --coverage-search --max-multihits 10 --b2-very-sensitive --min-coverage-intron ${min_intron_size} --max-coverage-intron ${max_intron_size} --microexon-search  --max-deletion-length 3 --max-insertion-length 3 --read-mismatches 2 --min-anchor $((min/3)) --splice-mismatches 0 --min-intron-length ${min_intron_size} --max-intron-length ${max_intron_size} --min-segment-intron ${min_intron_size} --max-segment-intron ${max_intron_size} --segment-mismatches=2 --segment-length=$((min/2)) --output-dir ${alignX_sample_out}/ --transcriptome-index ${tophat_transcriptome_index} -G ${refgff}"
		
		if [ -s "${alignL_sample_out}/all_junctions.junc" ]
		then   
			cmdtophat2="${cmdtophat2} --raw-juncs ${alignL_sample_out}/all_junctions.junc "
		fi
		
		cmdtophat2="${cmdtophat2} ${refbt2idx} ${fastqin}"
		
		if [ -e ${fastq2in} ] 
		then
			cmdtophat2="${cmdtophat2} ${fastq2in}"
		fi
		
		cmdtophat2="${cmdtophat2} &> ${alignX_sample_out}/${samplename}.tophat.log.txt"
		
		eval ${cmdtophat2}
		
		if [ -e "${alignX_sample_out}/junctions.bed" ]
		then
			bed_to_juncs < ${alignX_sample_out}/junctions.bed > ${alignX_sample_out}/junctions.junc 2> /dev/null
		else
			touch ${alignX_sample_out}/junctions.junc
		fi
		
		cat ${alignX_sample_out}/junctions.junc | nsort -u > ${alignX_sample_out}/all_junctions.junc
		
		samtools view -@ ${num_threads} -F 4 ${alignX_sample_out}/accepted_hits.bam | cut -f 1 | nsort -u > ${alignX_sample_out}/${samplename}.mapped.txt
		
		for fq1in in $(echo "${fastqin}" | tr "," "\n"); do
			rm -f ${alignX_sample_out}/${samplename}.unmapped_S*.fastq
			rm -f ${alignX_sample_out}/${samplename}.unmapped_U.fastq

			if [[ ${fq1in} =~ "prinseq_1_singletons" ]] || [[ ${fq1in} =~ "singleton_R1" ]]; then
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} >> ${alignX_sample_out}/${samplename}.unmapped_S1.fastq
			elif [[ ${fq1in} =~ "prinseq_2_singletons" ]] || [[ ${fq1in} =~ "singleton_R2" ]]; then
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} >> ${alignX_sample_out}/${samplename}.unmapped_S2.fastq
			else
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} >> ${alignX_sample_out}/${samplename}.unmapped_U.fastq
			fi
		
		done
		
		if [ -e ${fastq2in} ]
		then
			for fq2in in $(echo "${fastq2in}" | tr "," "\n"); do 
				if [[ ${fq2in} =~ "prinseq_1_singletons" ]] || [[ ${fq2in} =~ "singleton_R1" ]]; then
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_S1.fastq
				elif [[ ${fq2in} =~ "prinseq_2_singletons" ]] || [[ ${fq2in} =~ "singleton_R2" ]]; then
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_S2.fastq
				else
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_R2.fastq
				fi
			done
		fi
		
		#rm -fr ${alignX_sample_out}/tmp/
	}
	
	function bbmap_alignment_step {
		samplename=$1
		step=$2
		alignX_sample_out=$3
		alignL_sample_out=$4
		fastqin=$5
		fastq2in=$6
		
		bbmap_in2=""
		unmapped_second_end_fastq=""
		# (4) read unmapped
		samtools_accepted=" -F 4 " 
		samtools_unmapped=" -f 4 "
		
		if [ -e ${fastq2in} ] 
		then
			bbmap_fastqin2="in2=${fastq2in}"
			c=0
			# os singletons encontrados devem conter os correspondentes termos "null"
			for i in `echo ${fastqin} | tr "," "\n"`; do 
				if [ ${c} -gt 0 ]; then
					bbmap_fastqin2="${bbmap_fastqin2},null"
				fi
				c=$((c+1))
			done
			
			# (2) read mapped in proper pair 
			samtools_accepted=" -f 2 "
			samtools_unmapped=" -F 2 "
		fi
		
		cmdbbmap="bbwrap.sh -Xmx22g path=${refbbmapidx} build=1 in=${fastqin} ${bbmap_fastqin2} out=${alignX_sample_out}/all_hits.sam ambiguous=all trimq=0 mintrimlength=100 vslow=t strictmaxindel=t maxindel=3 minid=0.95 k=9 xstag=t mdtag=t nhtag=t xmtag=t amtag=t nmtag=t stoptag=t lengthtag=t idtag=t inserttag=t scoretag=t timetag=t boundstag=t xs=us sam=1.3 &> ${alignX_sample_out}/${samplename}.bbmap.log append=t"
		
		eval ${cmdbbmap}
		
		if [ ${step} == 3 ]
		then
			# (4) read unmapped
			samtools_accepted=" -F 4 "
			samtools_unmapped=" -f 4 "
		fi
		
		samtools view -@ ${num_threads} -b -S ${samtools_accepted} ${alignX_sample_out}/all_hits.sam > ${alignX_sample_out}/accepted_hits.bam 2> /dev/null
		samtools view -@ ${num_threads} -b -S ${samtools_unmapped} ${alignX_sample_out}/all_hits.sam > ${alignX_sample_out}/unmapped.bam 2> /dev/null
		
		rm -f ${alignX_sample_out}/all_hits.sam 
		
		samtools view -@ ${num_threads} -F 4 ${alignX_sample_out}/accepted_hits.bam | cut -f 1 | nsort -u > ${alignX_sample_out}/${samplename}.mapped.txt
		
		for fq1in in $(echo "${fastqin}" | tr "," "\n"); do
			if [[ ${fq1in} =~ "prinseq_1_singletons" ]]; then
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} > ${alignX_sample_out}/${samplename}.unmapped_S1.fastq
			elif [[ ${fq1in} =~ "prinseq_2_singletons" ]]; then
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} > ${alignX_sample_out}/${samplename}.unmapped_S2.fastq
			else
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} > ${alignX_sample_out}/${samplename}.unmapped_R1.fastq
			fi
		
		done
		
		if [ -e ${fastq2in} ]
		then
			for fq2in in $(echo "${fastq2in}" | tr "," "\n"); do 
				if [[ ${fq2in} =~ "prinseq_1_singletons" ]]; then
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_S1.fastq
				elif [[ ${fq2in} =~ "prinseq_2_singletons" ]]; then
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_S2.fastq
				else
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_R2.fastq
				fi
			done
		fi
		
	}
	
	function star_alignment_step {
		samplename=$1
		step=$2
		alignX_sample_out=$3
		alignL_sample_out=$4
		fastqin=$5
		fastq2in=$6
		
		star_fastq2in=""
		unmapped_second_end_fastq=""
		# (4) read unmapped
		samtools_accepted=" -F 4 " 
		samtools_unmapped=" -f 4 "
		
		if [ -e "$fastq2in" ]
		then
			star_fastq2in="$fastq2in"
			# (2) read mapped in proper pair 
			samtools_accepted=" -f 2 "
			samtools_unmapped=" -F 2 "
		fi
		
		cmdstar="STAR --limitBAMsortRAM 20000000000 --runThreadN ${num_threads} --outFilterMultimapNmax 10 --outFilterMismatchNmax 2 --outFilterMismatchNoverLmax 0.1 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique  --alignSJoverhangMin 35 --alignEndsType EndToEnd --alignIntronMin 5000000 --outSAMstrandField intronMotif --outSAMattributes All --outSAMprimaryFlag AllBestScore --bamRemoveDuplicatesType UniqueIdentical --outSAMtype BAM Unsorted --sjdbOverhang 100 --sjdbGTFfile ${refgff} --genomeDir ${refstaridx} --readFilesIn ${fastqin} ${star_fastq2in} --outFileNamePrefix ${alignX_sample_out}/ --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5 --alignMatesGapMax 600 --outSAMunmapped Within --chimOutType WithinBAM --chimSegmentMin 50"
		
		if [ -s "${alignL_sample_out}/all_junctions.junc" ]
		then   
			cmdstar="${cmdstar} --sjdbFileChrStartEnd ${alignL_sample_out}/all_junctions.junc "
		fi
		
		cmdstar="${cmdstar} &> ${alignX_sample_out}/${samplename}.STAR.log.txt"
		
		eval ${cmdstar}
		
		if [ -e "${alignX_sample_out}/SJ.out.tab" ]
		then
			perl -F"\t" -lane ' if (  ($F[4]>0) && ($F[5]>=0) && ($F[6]>=5) && ($F[7]>=0) && ($F[8]>10)  ) { print join("\t", @F[0,1,2],(($F[3]==1) ? "+" : (($F[3]==2) ? "-" : ".")) ); }' ${alignX_sample_out}/SJ.out.tab > ${alignX_sample_out}/junctions.junc
		else
			touch ${alignX_sample_out}/junctions.junc
		fi
		
		if [ -e "${alignL_sample_out}/all_junctions.junc" ]; then
			cat ${alignL_sample_out}/all_junctions.junc ${alignX_sample_out}/junctions.junc | nsort -u > ${alignX_sample_out}/all_junctions.junc
		else
			cat ${alignX_sample_out}/junctions.junc | nsort -u > ${alignX_sample_out}/all_junctions.junc
		fi
		
		if [ ${step} == 3 ]
		then
			# (4) read unmapped
			samtools_accepted=" -F 4 "
			samtools_unmapped=" -f 4 "
		fi
		
		samtools view -@ ${num_threads} -b ${samtools_accepted} ${alignX_sample_out}/Aligned.out.bam > ${alignX_sample_out}/accepted_hits.bam 2> /dev/null
		samtools view -@ ${num_threads} -h ${samtools_unmapped} ${alignX_sample_out}/Aligned.out.bam > ${alignX_sample_out}/unmapped.bam 2> /dev/null
		samtools view -@ ${num_threads} ${alignX_sample_out}/accepted_hits.bam | cut -f 1 | nsort -u > ${alignX_sample_out}/${samplename}.mapped.txt
		
		for fq1in in $(echo "${fastqin}" | tr "," "\n"); do
			if [[ ${fq1in} =~ "prinseq_1_singletons" ]]; then
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} > ${alignX_sample_out}/${samplename}.unmapped_S1.fastq
			elif [[ ${fq1in} =~ "prinseq_2_singletons" ]]; then
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} > ${alignX_sample_out}/${samplename}.unmapped_S2.fastq
			else
				cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq1in} > ${alignX_sample_out}/${samplename}.unmapped_R1.fastq
			fi
		
		done
		
		if [ -e ${fastq2in} ]
		then
			for fq2in in $(echo "${fastq2in}" | tr "," "\n"); do 
				if [[ ${fq2in} =~ "prinseq_1_singletons" ]]; then
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_S1.fastq
				elif [[ ${fq2in} =~ "prinseq_2_singletons" ]]; then
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_S2.fastq
				else
					cat ${alignX_sample_out}/${samplename}.mapped.txt | cut -f 1 | nsort -u | pullseq -e -N -i ${fq2in} >> ${alignX_sample_out}/${samplename}.unmapped_R2.fastq
				fi
			done
		fi
	
	}
	
	echo "* Pre-processing:"
	
	for fastq in `ls ${fastx_out}/*_R1.fastq*`; do
				
		# fqname - string com o nome do arquivo fastq
		fqlink=$(readlink -f ${fastq})
		fqtype=`file -b -i ${fqlink} | cut -d ';' -f 1`
		if [ ${fqtype} == "text/plain" ]; then
			fqname=`basename ${fastq} _R1.fastq`
		else
			fqname=`basename ${fastq} _R1.fastq.gz`
		fi
		
		mate_fastq=`echo ${fastq} | sed 's/_R1.fastq/_R2.fastq/'`
		
		echo "   ${fqname} ..."	
		
		# Execução do align2 para obter o alinhamento de cada réplica biológica em relação ao genoma
		# caso ainda não tenha sido executado
		
		align0_sample_out="${align0_out}/${fqname}"
		mkdir -p ${align0_sample_out}
		
		if [ !  -e "${align0_sample_out}/accepted_hits.bam" ]
		then
			if [ ! ${fqtype} == "text/plain" ]; then
				fastq=`echo ${fastq} | sed 's/_R1.fastq.gz$/_R1.fastq/'`;
				mate_fastq=`echo ${mate_fastq} | sed 's/_R2.fastq.gz$/_R2.fastq/'`;
				echo " Uncompressing ${fastq} ..."
				gunzip ${fastq}
				echo " Uncompressing ${mate_fastq}..."
				gunzip ${mate_fastq}
			fi
			
			echo "      Aligning ${fqname} against genome - step 0/4 (${aligner}) ..."
			
			if [ ${aligner} == "tophat2" ]; then
				
				tophat_alignment_step ${fqname} 0 ${align0_sample_out} ${align0_sample_out} ${fastq} ${mate_fastq}
				
			elif [ ${aligner} == "bbmap" ]; then
				
				bbmap_alignment_step ${fqname} 0 ${align0_sample_out} ${align0_sample_out} ${fastq} ${mate_fastq}

			elif [ ${aligner} == "star" ]; then

				star_alignment_step ${fqname} 0 ${align0_sample_out} ${align0_sample_out} ${fastq} ${mate_fastq}
			
			elif [ ${aligner} == "hisat2" ]; then

				hisat2_alignment_step ${fqname} 0 ${align0_sample_out} ${align0_sample_out} ${fastq} ${mate_fastq}
			else
				echo "Wrong aligner (${aligner})"
				exit
			fi
			
			if [ -e ${fastq} ]; then
				if [ ! -L ${fastq} ]; then
					echo "         Compressing ${fqname} [raw2] ..."
					gzip ${fastq} ${mate_fastq}
				fi
			fi
		fi
		
		if [ ! -e "${seqpurge_out}/${fqname}.seqpurge_R1.fastq" ] 
		then   
			mkdir -p ${seqpurge_out}
			
			if [ -e ${align0_sample_out}/${fqname}.unmapped_R1.fastq ] && [ ${align0_sample_out}/${fqname}.unmapped_R1.fastq ]
			then

				# seqpurge_adapter_r1
				if [ ! ${seqpurge_adapter_r1} ]; then
					echo "seqpurge_adapter_r1 not defined";
					exit;
				fi

				# seqpurge_adapter_r2
				if [ ! ${seqpurge_adapter_r2} ]; then
					echo "seqpurge_adapter_r2 not defined";
					exit;
				fi
				
				echo "        Adapter trimming - step 1/3 (SeqPurge) [R1 & R2] ..."
				
				SeqPurge \
					-in1 ${align0_sample_out}/${fqname}.unmapped_R1.fastq \
					-in2 ${align0_sample_out}/${fqname}.unmapped_R2.fastq \
					-out1 ${seqpurge_out}/${fqname}.seqpurge_R1.fastq.gz \
					-out2 ${seqpurge_out}/${fqname}.seqpurge_R2.fastq.gz \
					-a1 ${seqpurge_adapter_r1} \
					-a2 ${seqpurge_adapter_r2} \
					-summary ${seqpurge_out}/${fqname}.summary_file.txt \
					-qwin 3 \
					-min_len 25 \
					-out3 ${seqpurge_out}/${fqname}.seqpurge_singleton \
					-threads ${num_threads} \
					-ncut 3 \
					-progress
					
					gunzip ${seqpurge_out}/${fqname}.seqpurge_R1.fastq.gz
					gunzip ${seqpurge_out}/${fqname}.seqpurge_R2.fastq.gz
					
					if [ ${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq.gz ]; then
						gunzip ${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq.gz
					fi
					
					if [ ${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq.gz ]; then
						gunzip ${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq.gz
					fi
			else
				echo "        Adapter trimming - step 1/3 (SeqPurge) [R1 & R2] ... Not performed!"
				ln -s ${align0_sample_out}/${fqname}.unmapped_R1.fastq ${seqpurge_out}/${fqname}.seqpurge_R1.fastq
				ln -s ${align0_sample_out}/${fqname}.unmapped_R2.fastq ${seqpurge_out}/${fqname}.seqpurge_R2.fastq
			fi
		fi
		
		align1_sample_out="${align1_out}/${fqname}"
		mkdir -p ${align1_sample_out}

		if [ !  -e "${align1_sample_out}/accepted_hits.bam" ]
		then
			echo "      Aligning ${fqname} against genome - step 1/4 (${aligner}) ..."
			
			if [ ${aligner} == "tophat2" ]; then

				tophat_alignment_step ${fqname} 1 ${align1_sample_out} ${align0_sample_out} "${seqpurge_out}/${fqname}.seqpurge_R1.fastq" "${seqpurge_out}/${fqname}.seqpurge_R2.fastq"

			elif [ ${aligner} == "bbmap" ]; then
				
				bbmap_alignment_step ${fqname} 1 ${align1_sample_out} ${align0_sample_out} "${seqpurge_out}/${fqname}.seqpurge_R1.fastq" "${seqpurge_out}/${fqname}.seqpurge_R2.fastq"

			elif [ ${aligner} == "star" ]; then
		
				star_alignment_step ${fqname} 1 ${align1_sample_out} ${align0_sample_out} "${seqpurge_out}/${fqname}.seqpurge_R1.fastq" "${seqpurge_out}/${fqname}.seqpurge_R2.fastq"
			elif [ ${aligner} == "hisat2" ]; then
		
				hisat2_alignment_step ${fqname} 1 ${align1_sample_out} ${align0_sample_out} "${seqpurge_out}/${fqname}.seqpurge_R1.fastq" "${seqpurge_out}/${fqname}.seqpurge_R2.fastq"
			else 
				echo "Wrong aligner (${aligner})"
				exit
			fi
			
		fi
		
		if [ ! -e "${scythe_out}/${fqname}.scythe_R1.fastq" ] && [ ! -e "${scythe_out}/${fqname}.scythe_R1.fastq.gz" ]
		then   

			mkdir -p ${scythe_out}
			
			echo "         Adapter trimming - step 2/3 (scythe) [R1] ..."
			
			# -m ${scythe_out}/${fqname}.scythe.matches_R1.txt
			scythe -p 0.1 -o ${scythe_out}/${fqname}.scythe_R1.fastq -n 5 -M 0  -q sanger  -a ${scythe_adapters_r1_path} -i ${align1_sample_out}/${fqname}.unmapped_R1.fastq 1> ${scythe_out}/${fqname}.scythe_R1.log.txt 2> ${scythe_out}/${fqname}.scythe_R1.err.txt

			cat ${scythe_out}/${fqname}.scythe_R1.fastq | sed 's/ *$//' > ${scythe_out}/${fqname}.scythe_R1.tmp
			
			mv ${scythe_out}/${fqname}.scythe_R1.tmp ${scythe_out}/${fqname}.scythe_R1.fastq
		fi
		
		if [ -e "${seqpurge_out}/${fqname}.seqpurge_R2.fastq" ]
		then
			if [ ! -e "${scythe_out}/${fqname}.scythe_R2.fastq" ] && [ ! -e "${scythe_out}/${fqname}.scythe_R2.fastq.gz" ]
			then   

				echo "         Adapter trimming - step 2/3 (scythe) [R2] ..."
				
				# -m ${scythe_out}/${fqname}.scythe.matches_R2.txt
				scythe -p 0.1 -o ${scythe_out}/${fqname}.scythe_R2.fastq -n 5 -M 0  -q sanger  -a ${scythe_adapters_r2_path} -i ${align1_sample_out}/${fqname}.unmapped_R2.fastq 1> ${scythe_out}/${fqname}.scythe_R2.log.txt 2> ${scythe_out}/${fqname}.scythe_R2.err.txt

				cat ${scythe_out}/${fqname}.scythe_R2.fastq | sed 's/ *$//' > ${scythe_out}/${fqname}.scythe_R2.tmp
			
				mv ${scythe_out}/${fqname}.scythe_R2.tmp ${scythe_out}/${fqname}.scythe_R2.fastq
			fi
		fi

		# Execução do align2 para obter o alinhamento de cada réplica biológica em relação ao genoma
		# caso ainda não tenha sido executado
		
		align2_sample_out="${align2_out}/${fqname}"

		mkdir -p ${align2_sample_out}

		if [ !  -e "${align2_sample_out}/accepted_hits.bam" ]
		then
			echo "      Aligning ${fqname} against genome - step 2/4 (${aligner}) ..."
			
			if [ ${aligner} == "tophat2" ]; then

				tophat_alignment_step ${fqname} 2 ${align2_sample_out} ${align1_sample_out} "${scythe_out}/${fqname}.scythe_R1.fastq" "${scythe_out}/${fqname}.scythe_R2.fastq"

			elif [ ${aligner} == "bbmap" ]; then
				
				bbmap_alignment_step ${fqname} 2 ${align2_sample_out} ${align1_sample_out} "${scythe_out}/${fqname}.scythe_R1.fastq" "${scythe_out}/${fqname}.scythe_R2.fastq"

			elif [ ${aligner} == "star" ]; then
		
				star_alignment_step ${fqname} 2 ${align2_sample_out} ${align1_sample_out} "${scythe_out}/${fqname}.scythe_R1.fastq" "${scythe_out}/${fqname}.scythe_R2.fastq"
			elif [ ${aligner} == "hisat2" ]; then
		
				hisat2_alignment_step ${fqname} 2 ${align2_sample_out} ${align1_sample_out} "${scythe_out}/${fqname}.scythe_R1.fastq" "${scythe_out}/${fqname}.scythe_R2.fastq"	
			else 
				echo "Wrong aligner (${aligner})"
				exit
			fi
			
			if [ -e "${scythe_out}/${fqname}.scythe_R1.fastq"  ]; then
				echo "         Compressing ${fqname} [scythe] ..."
				gzip ${scythe_out}/${fqname}.*
			fi
		fi

		if [ ! -e "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" ] && [ ! -e "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq.gz" ]
		then
			
			mkdir -p ${cutadapt_out}

			echo "         Adapter trimming - step 3/3 (cutadapt) [R1] ..."
			
			cutadapt_5p_adapter_r1=`perl -ane 'chomp; if ($_=~/^>(\S+)/) { } else { print "--front=\"$_\""," "; }' ${cutadapt_5p_adapters_r1_path}`;

			cutadapt --format=fastq ${cutadapt_5p_adapter_r1} --error-rate=0.1 --times=2 --overlap=5 --minimum-length=0 --output=${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq  ${align2_sample_out}/${fqname}.unmapped_R1.fastq > ${cutadapt_out}/${fqname}.scythe.cutadapt_R1.report.txt 2> ${cutadapt_out}/${fqname}.scythe.cutadapt_R1.log.txt

		fi
		if [ -e "${align2_sample_out}/${fqname}.unmapped_R2.fastq" ] || [ -e "${align2_sample_out}/${fqname}.unmapped_R2.fastq.gz" ] 
		then
			if [ ! -e "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq" ] && [ ! -e "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq.gz" ]
			then
				echo "         Adapter trimming - step 3/3 (cutadapt) [R2] ..."
			
				cutadapt_5p_adapter_r2=`perl -ane 'chomp; if ($_=~/^>(\S+)/) { } else { print "--front=\"$_\""," "; }' ${cutadapt_5p_adapters_r2_path}`;

				cutadapt --format=fastq ${cutadapt_5p_adapter_r2} --error-rate=0.1 --times=2 --overlap=5 --minimum-length=0 --output=${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq  ${align2_sample_out}/${fqname}.unmapped_R2.fastq > ${cutadapt_out}/${fqname}.scythe.cutadapt_R2.report.txt 2> ${cutadapt_out}/${fqname}.scythe.cutadapt_R2.log.txt

			fi
		fi
		
		
		# Execução do align3 para obter o alinhamento de cada réplica biológica em relação ao genoma
		# caso ainda não tenha sido executado
		
		align3_sample_out="${align3_out}/${fqname}"
		mkdir -p ${align3_sample_out}
		
		if [ !  -e "${align3_sample_out}/accepted_hits.bam" ]
		then
			echo "      Aligning ${fqname} against genome - step 3/4 (${aligner}) ..."
			
			if [ ${aligner} == "tophat2" ]; then
				
				tophat_alignment_step ${fqname} 3 ${align3_sample_out} ${align2_sample_out} "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq"
				
			elif [ ${aligner} == "bbmap" ]; then
				
				bbmap_alignment_step ${fqname} 3 ${align3_sample_out} ${align2_sample_out} "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq"
				
			elif [ ${aligner} == "star" ]; then
				
				star_alignment_step ${fqname} 3 ${align3_sample_out} ${align2_sample_out} "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq"
			
			elif [ ${aligner} == "hisat2" ]; then
				
				hisat2_alignment_step ${fqname} 3 ${align3_sample_out} ${align2_sample_out} "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq" "${cutadapt_out}/${fqname}.scythe.cutadapt_R2.fastq"
		
			else 
				echo "Not found aligner (${aligner})"
				exit
			fi
			
			if [ -e "${cutadapt_out}/${fqname}.scythe.cutadapt_R1.fastq"  ]; then
				echo "         Compressing ${fqname} [cutadapt] ..."
				gzip ${cutadapt_out}/${fqname}.*
			fi
			
		fi

		if [ ! -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ] && [ ! -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ] && [ ! -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq.gz" ] && [ ! -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq.gz" ]
		then

			echo "         Quality and poly-A/T trimming (prinseq) ..."

			mkdir -p ${prinseq_out}
			

			fastq2in=""
			
			if [ -e "${align3_sample_out}/${fqname}.unmapped_R2.fastq" ]
			then
				fastq2in="-fastq2 ${align3_sample_out}/${fqname}.unmapped_R2.fastq"
			fi

			# trimming reads by quality value
			prinseq-lite.pl -fastq ${align3_sample_out}/${fqname}.unmapped_R1.fastq ${fastq2in} -out_format 3 -out_good ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq -out_bad ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq-bad -qual_noscale -no_qual_header -min_len 15 -ns_max_p 80 -noniupac -trim_qual_right 25 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 3 -trim_qual_step 1 &> ${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.log.txt

		fi
		
		# Execução do align4 para obter o alinhamento de cada réplica biológica em relação ao genoma
		# caso ainda não tenha sido executado
		
		align4_sample_out="${align4_out}/${fqname}"
		mkdir -p ${align4_sample_out}

		if [ !  -e "${align4_sample_out}/accepted_hits.bam" ]
		then
			echo "      Aligning ${fqname} against genome - step 4/4 (${aligner}) ..."

			if [ ${aligner} == "tophat2" ]
			then
				
				tophat_fastq_in=""
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]
				then
					tophat_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq"
				elif [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ]
				then
					tophat_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq"
					if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ] &&
					   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ]
					then
						tophat_fastq_in="${tophat_fastq_in},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq"
					fi
					if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ] &&
					   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ] 
					then
						tophat_fastq_in="${tophat_fastq_in},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq"
					fi
				else
					echo "Missing input file for tophat2 alignment step 3"
				fi

				if [ -e "${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq" ] &&
				   [ -s "${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq" ] 
				then
					tophat_fastq_in="${tophat_fastq_in},${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq"
				fi
				if [ -e "${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq" ] &&
				   [ -s "${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq" ]
				then
					tophat_fastq_in="${tophat_fastq_in},${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq"
				fi
				
				tophat_alignment_step ${fqname} 4 ${align4_sample_out} ${align3_sample_out} ${tophat_fastq_in} "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq"

			elif [ ${aligner} == "hisat2" ]
			then
				
				hisat2_fastq_in=""
				hisat2_fastq_u=""
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq.gz" ] || 
				   [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq.gz" ]; then
					gunzip ${prinseq_out}/${fqname}.*.fastq.gz
					gunzip ${prinseq_out}/${fqname}.*.log.txt.gz
				fi

				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ] &&
				   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]
				then
					hisat2_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq"
				elif [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ]
				then
					hisat2_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq"
					if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ] &&
					   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ]
					then
						hisat2_fastq_u="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq"
					fi
					if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ] &&
					   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ]
					then
						if [ ${hisat2_fastq_u} ]; then
							hisat2_fastq_u="${hisat2_fastq_u},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq"
						else
							hisat2_fastq_u="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq"
						fi
					fi
				else
					echo "Missing input file for hisat2 alignment step 4"
				fi
				
				if [ -e "${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq" ] &&
				   [ -s "${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq" ]
				then
					if [ ${hisat2_fastq_u} ]; then
						hisat2_fastq_u="${hisat2_fastq_u},${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq"
					else
						hisat2_fastq_u="${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq"
					fi
				fi
				if [ -e "${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq" ] &&
				   [ -s "${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq" ]
				then
					if [ ${hisat2_fastq_u} ]; then
						hisat2_fastq_u="${hisat2_fastq_u},${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq"
					else
						hisat2_fastq_u="${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq"
					fi
				fi

				hisat2_alignment_step ${fqname} 4 ${align4_sample_out} ${align3_sample_out} ${hisat2_fastq_in} "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq" ${hisat2_fastq_u}

			elif [ ${aligner} == "bbmap" ]; then
				
				bbmap_fastq_in=""
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ] &&
				   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]
				then
					bbmap_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq"
				elif [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ] &&
				     [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ]
				then
					bbmap_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq"
					if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ]; then
						bbmap_fastq_in="${bbmap_fastq_in},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq"
					fi
					if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ]; then
						bbmap_fastq_in="${bbmap_fastq_in},${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq"
					fi
				else
					echo "Missing input file for bbmap alignment step 3"
				fi
				
				if [ -e "${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq" ] &&
				   [ -s "${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq" ]
				then
					bbmap_fastq_in="${bbmap_fastq_in},${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq"
				fi
				if [ -e "${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq" ] &&
				   [ -s "${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq" ]
				then
					bbmap_fastq_in="${bbmap_fastq_in},${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq"
				fi
				
				bbmap_alignment_step ${fqname} 4 ${align4_sample_out} ${align3_sample_out} ${bbmap_fastq_in} "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq"
				
			else
				star_fastq_in=""
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ] &&
				   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]
				then
					star_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq"
				elif [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ] &&
				     [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ]
				then
					star_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq"
				else
					echo "Missing input file for star alignment step 3"
				fi
				
				star_alignment_step ${fqname} 4 ${align4_sample_out} ${align3_sample_out} ${star_fastq_in} "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2.fastq"
				
				### seqpurge singletons
				if [ -e "${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq" ] &&
				   [ -s "${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq" ]; then
					star_fastq_in="${seqpurge_out}/${fqname}.seqpurge_singleton_R1.fastq"
					mkdir -p ${align4_sample_out}/seqpurge_singleton1
					star_alignment_step ${fqname} 4 ${align4_sample_out}/seqpurge_singleton1 ${align4_sample_out} ${star_fastq_in} ""
					
				fi
				if [ -e "${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq" ] &&
				   [ -s "${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq" ]; then
					star_fastq_in="${seqpurge_out}/${fqname}.seqpurge_singleton_R2.fastq"
					mkdir -p ${align4_sample_out}/seqpurge_singleton2
					star_alignment_step ${fqname} 4 ${align4_sample_out}/seqpurge_singleton2 ${align4_sample_out} ${star_fastq_in} ""
				fi
			
				if [ ! -e "${align4_sample_out}/header.sam" ]
				then 
					samtools view -@ ${num_threads} -H ${align4_sample_out}/accepted_hits.bam > ${align4_sample_out}/header.sam 2> /dev/null
				fi
				
				if [ -e "${align4_sample_out}/seqpurge_singleton1/accepted_hits.bam" ]; then
					mv ${align4_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.tmp

					samtools merge -@ ${num_threads} -h ${align4_sample_out}/header.sam ${align4_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.tmp ${align4_sample_out}/seqpurge_singleton1/accepted_hits.bam 2> /dev/null
					
					rm -f ${align4_sample_out}/accepted_hits.tmp
				fi
				
				if [ -e "${align4_sample_out}/seqpurge_singleton2/accepted_hits.bam" ]; then
					mv ${align4_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.tmp

					samtools merge -@ ${num_threads} -h ${align4_sample_out}/header.sam ${align4_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.tmp ${align4_sample_out}/seqpurge_singleton2/accepted_hits.bam 2> /dev/null
					
					rm -f ${align4_sample_out}/accepted_hits.tmp
				fi
				
				### prinseq singletons
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ] && 
				   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq" ]; then
					star_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1_singletons.fastq"
					mkdir -p ${align4_sample_out}/singleton1
					star_alignment_step ${fqname} 4 ${align4_sample_out}/singleton1 ${align4_sample_out} ${star_fastq_in} ""
					
				fi
				if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ] && 
				   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq" ]; then
					star_fastq_in="${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_2_singletons.fastq"
					mkdir -p ${align4_sample_out}/singleton2
					star_alignment_step ${fqname} 4 ${align4_sample_out}/singleton2 ${align4_sample_out} ${star_fastq_in} ""
				fi
			
				if [ ! -e "${align4_sample_out}/header.sam" ]
				then 
					samtools view -@ ${num_threads} -H ${align4_sample_out}/accepted_hits.bam > ${align4_sample_out}/header.sam 2> /dev/null
				fi
				
				if [ -e "${align4_sample_out}/singleton1/accepted_hits.bam" ]; then
					mv ${align4_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.tmp

					samtools merge -@ ${num_threads} -h ${align4_sample_out}/header.sam ${align4_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.tmp ${align4_sample_out}/singleton1/accepted_hits.bam 2> /dev/null
					
					rm -f ${align4_sample_out}/accepted_hits.tmp
				fi
				
				if [ -e "${align4_sample_out}/singleton2/accepted_hits.bam" ]; then
					mv ${align4_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.tmp

					samtools merge -@ ${num_threads} -h ${align4_sample_out}/header.sam ${align4_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.tmp ${align4_sample_out}/singleton2/accepted_hits.bam 2> /dev/null
					
					rm -f ${align4_sample_out}/accepted_hits.tmp
				fi

			fi
			
			if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ] && 
			   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq.fastq" ]; then
				echo "         Compressing ${fqname} [prinseq] ..."
				gzip ${prinseq_out}/${fqname}.*	
			fi
			if [ -e "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ] &&
			   [ -s "${prinseq_out}/${fqname}.scythe.cutadapt.prinseq_1.fastq" ]; then
				echo "         Compressing ${fqname} [prinseq] ..."
				gzip ${prinseq_out}/${fqname}.*
			fi
			
		fi

		finalalign_sample_out="${finalalign_out}/${fqname}"
		mkdir -p ${finalalign_sample_out}

		if [ ! -e "${finalalign_sample_out}/accepted_hits.bam" ]
		then
			echo "      Merging alignments of ${fqname} against genome (samtools) ..."

			if [ ! -e "${finalalign_sample_out}/header.sam" ]
			then 
				samtools view -@ ${num_threads} -H ${align0_sample_out}/accepted_hits.bam > ${finalalign_sample_out}/header.sam 2> /dev/null
			fi

			samtools merge -@ ${num_threads} -h ${finalalign_sample_out}/header.sam ${finalalign_sample_out}/accepted_hits.bam ${align0_sample_out}/accepted_hits.bam ${align1_sample_out}/accepted_hits.bam ${align2_sample_out}/accepted_hits.bam ${align3_sample_out}/accepted_hits.bam ${align4_sample_out}/accepted_hits.bam 2> /dev/null

			samtools sort -@ ${num_threads} ${finalalign_sample_out}/accepted_hits.bam -o ${finalalign_sample_out}/accepted_hits_sorted.bam

			mv ${finalalign_sample_out}/accepted_hits_sorted.bam ${finalalign_sample_out}/accepted_hits.bam

		fi
		
		
		
		
	done
fi

if [ ${execopt} == "diff" ] || [ ${execopt} == "all" ]; then
	
	
	for K in "${!MYDESIGN[@]}"; do
		

        	echo "Processing design: $K ...";

		if [ ! -e ${MYDESIGN[$K]} ]; then
			echo "Missing design file ${MYDESIGN[$K]}"
			exit
		fi
		
		if [ ! -d ${fastx_out} ]; then
			echo "Error: Not found $fastx_out Please, you need to run align first!"
			exit
		fi
		
		echo "* Creating directories (differential analysis step)"
		
		# cufflinks_out - caminho para o diretório de saída do cufflinks
		cufflinks_out="${basedir_out}/$K/cufflinks"
		
		# cuffmerge_out - caminho para o diretório de saída do cuffmerge
		cuffmerge_out="${basedir_out}/$K/cuffmerge"
		
		# cuffquant_out - caminho para o diretório de saída do cuffquant
		cuffquant_out="${basedir_out}/$K/cuffquant"
		
		# cuffnorm_out - caminho para o diretório de saída do cuffnorm
		cuffnorm_out="${basedir_out}/$K/cuffnorm"
		
		# cuffdiff_out - caminho para o diretório de saída do cuffdiff
		cuffdiff_out="${basedir_out}/$K/cuffdiff"
		# Execução do cuffmerge para fusionar os transcriptomas gerados em um transcriptoma referência
		
		mkdir -p ${cuffmerge_out}
		mkdir -p ${cuffquant_out}
		mkdir -p ${cuffnorm_out}
		mkdir -p ${cuffdiff_out}
		
		if [ ${cuffdiff_direct} == "yes" ]; then
			
			if [ ! -e ${cuffmerge_out}/merged.gtf ]; then
				ln -s ${refgff} ${cuffmerge_out}/merged.gtf
			fi
			
		else
			if [ ! ${cuffdiff_direct} == "no" ]; then
				echo "Error found! Parameter value for cufflinks_direct not recognized (${cuffdiff_direct})"
				exit;
			fi
		fi

		if [ ! -e "${cuffmerge_out}/merged.gtf" ]
		then
			
			mkdir -p ${cufflinks_out}
			
			if [ ${cufflinks_guided} == "yes" ]; then
				cufflinks_gtf_option="--GTF-guide"
			else
				if [ ${cufflinks_guided} == "no" ]; then
					cufflinks_gtf_option="--GTF"
				else
					echo "Error found! Parameter value for cufflinks_guided not recognized (${cufflinks_guied})"
					exit;
				fi
			fi
		
			for fastq in `ls ${fastx_out}/*_R1.fastq*`; do
				
				# fqname - string com o nome do arquivo fastq
				fqlink=$(readlink -f ${fastq})
				fqtype=`file -b -i ${fqlink} | cut -d ';' -f 1`
				if [ ${fqtype} == "text/plain" ]; then
					fqname=`basename ${fastq} _R1.fastq`
				else
					fqname=`basename ${fastq} _R1.fastq.gz`
				fi
				
				mate_fastq=`echo ${fastq} | sed 's/_R1.fastq.gz$/_R2.fastq.gz/'`
				
				echo "   ${fqname} ..."	
				# Execução do cufflinks para obter a montagem para cada réplica biológica
				
				finalalign_sample_out="${finalalign_out}/${fqname}"
				if [ ! -e ${finalalign_sample_out} ]; then
					echo "Directory ${finalalign_sample_out} doesn't exist"
					exit
				fi
				# caso ainda não tenha sido executado
				cufflinks_sample_out="${cufflinks_out}/${fqname}"
				if [ !  -e "${cufflinks_sample_out}/transcripts.gtf" ]
				then
						
					cufflinks_gtf_cmdline="${cufflinks_gtf_option} ${refgff}"
					if [ ${cufflinks_gtf} ]; then
						if [ -e ${cufflinks_gtf} ]; then
							cufflinks_gtf_cmdline="${cufflinks_gtf_option} ${cufflinks_gtf}"
						else
							if [ ${cufflinks_gtf} = "no" ]; then
								cufflinks_gtf_cmdline=""		
							else
								echo "Wrong parameter (cufflinks_gtf) value (${cufflinks_gtf})"
								exit
							fi
						fi				
					fi
					mkdir -p ${cufflinks_sample_out}
									
					echo "      Assembling ${fqname} (cufflinks) ..."
					
					cufflinks --num-threads ${num_threads} --library-type ${cufflinks_library_type} --multi-read-correct ${cufflinks_gtf_cmdline} --frag-bias-correct ${refgenomefa} --total-hits-norm --output-dir ${cufflinks_sample_out} ${finalalign_sample_out}/accepted_hits.bam --max-bundle-frags 2000000 &> ${cufflinks_sample_out}/${fqname}.cufflinks.log.txt
				fi
				
			done
			
			# removendo arquivo com o caminho de cada transcriptoma (.gtf)
			rm -f ${cuffmerge_out}/assembly_GTF_list.txt
			
			for bs in $(cut -f 1 ${MYDESIGN[$K]}); do
				# caso ainda não tenha sido executado
				cufflinks_sample_out="${cufflinks_out}/${bs}"
				if [ ! -d ${cufflinks_sample_out} ]; then
					echo "Wrong design $K in file ${MYDESIGN[$K]}"
					rm -f ${cuffmerge_out}/assembly_GTF_list.txt
					exit
				else
					# gravando no arquivo ${cuffmerge_out}/assembly_GTF_list.txt o caminho de cada transcriptoma para cada réplica biológica. Esse arquivo é requerido para a execução do cuffmerge
					echo "${cufflinks_sample_out}/transcripts.gtf" >> ${cuffmerge_out}/assembly_GTF_list.txt
				fi
			done

			cuffmerge_gtf_cmdline="--ref-gtf ${refgff}"
			if [ ${cuffmerge_gtf} ]; then
				if [ -e ${cuffmerge_gtf} ]; then
					cuffmerge_gtf_cmdline="--ref-gtf ${cuffmerge_gtf}"
				else
					if [ ${cuffmerge_gtf} = "no" ]; then
						cuffmerge_gtf_cmdline=""		
					else
						echo "Wrong parameter (cuffmerge_gtf) value (${cuffmerge_gtf})"
						exit
					fi
				fi				
			fi
			
			rowcount=`wc -l ${cuffmerge_out}/assembly_GTF_list.txt | sed 's/ .*//'`;
			echo "* Merging ${rowcount} transcriptomes (cuffmerge) ..."

			cuffmerge --min-isoform-fraction 0.05 --num-threads ${num_threads} --ref-sequence ${refgenomefa} -o ${cuffmerge_out} ${cuffmerge_gtf_cmdline} ${cuffmerge_out}/assembly_GTF_list.txt &> ${cuffmerge_out}/cuffmerge.log.txt

		fi
		
		echo "Creating quantification files (.cxb):"
		for bs in $(cut -f 1 ${MYDESIGN[$K]}); do
			
			echo "	${bs} ..."
			
			mkdir -p ${cuffquant_out}/${bs}
			
			if [ ! -e "${cuffquant_out}/${bs}/abundances.cxb" ]
			then
				cuffquant --num-threads ${num_threads} --library-type ${cufflinks_library_type} --multi-read-correct --frag-bias-correct ${refgenomefa} ${cuffmerge_out}/merged.gtf ${finalalign_out}/${bs}/accepted_hits.bam -o ${cuffquant_out}/${bs} --max-bundle-frags 999999999 &> ${cuffquant_out}/${bs}.log.txt
			fi
			
		done;
		
		biogroupscxb=()
		biogroupslabel=($(cut -f 2 ${MYDESIGN[$K]} | sort -u | xargs) )
		echo "Checking cxb files:"
		for bg in ${biogroupslabel[@]}; do
			
			echo "	${bg}"
				
			cxblist=()
			
			for br in $( grep "${bg}" ${MYDESIGN[$K]} | cut -f 1); do
				echo "		${br}"
				if [ ! -e ${cuffquant_out}/${br}/abundances.cxb ]; then
					echo "Error: not found ${cuffquant_out}/${br}/abundances.cxb"
					exit
				fi
				cxblist=( ${cxblist[@]} ${cuffquant_out}/${br}/abundances.cxb )
			done
			
			# aqui temos que usar o asterisco ("*") para expandir o conteúdo de cxblist pois queremos apenas uma palavra, separada pela variável especial IFS e não múltiplas palavras como quando usamos o arroba ("@")
			biogroupscxb=(${biogroupscxb[@]} $(IFS=, ; echo "${cxblist[*]}"))
		done;
		
		if [ ! -e "${cuffnorm_out}/run.info" ]
		then
			echo "* Scaling transcriptomes (cuffnorm) ..."
			
			cuffnorm --num-threads ${num_threads} --library-type ${cufflinks_library_type} --library-norm-method geometric --output-format simple-table --compatible-hits-norm --output-dir ${cuffnorm_out} --labels $(IFS=, ; echo "${biogroupslabel[*]}") ${cuffmerge_out}/merged.gtf ${biogroupscxb[@]} &> ${cuffnorm_out}/cuffnorm.log

		fi 

		# Execução do cuffdiff para obter as estimativas de expressão gênica diferencial
		# caso não tenha sido executado
		if [ ! -e "${cuffdiff_out}/run.info" ]
		then
			echo "* Estimating Differential Expression (cuffdiff) among groups:" "${biogroupslabel[@]}" 
			cuffdiff --num-threads ${num_threads} --library-type ${cufflinks_library_type} --multi-read-correct --frag-bias-correct ${refgenomefa} --library-norm-method geometric --compatible-hits-norm --output-dir ${cuffdiff_out} ${cuffmerge_out}/merged.gtf ${biogroupscxb[@]} --FDR  0.05 --labels $(IFS=, ; echo "${biogroupslabel[*]}") --max-bundle-frags 999999999 &> ${cuffdiff_out}/cuffdiff.log
		fi

		
	done
		

fi
	
# Documentação do bash
#	http://www.gnu.org/software/bash/manual/bashref.html#Arrays
# Discussão sobre * ou @ ao expandir valores do array
#	Suppose we have perls=(perl-one perl-two) and use ls "${perls[*]}" -- that's equivalent to ls "perl-one perl-two", which will look for single file named perl-one perl-two, which is probably not what you wanted. ls "${perls[@]}" is equivalent to ls "perl-one" "perl-two", which is much more likely to do something useful.
# 	http://www.coderexception.com/0166bbz36XXQJQQX/a-confusion-about-array-versus-array-in-the-context-of-a-bash-completion
