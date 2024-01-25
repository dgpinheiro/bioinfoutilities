#!/bin/bash

# ALINHADOR tophat OU star
aligner=$1

if [ ! ${aligner} ]; then
	echo "[ERROR] Missing aligner (tophat or star)." 1>&2
	exit
fi

if [ "${aligner}" != "tophat" ] &&
   [ "${aligner}" != "star" ]; then

	echo "[ERROR] Aligner must be \"tophat\" or \"star\" (${aligner})." 1>&2
  	exit
fi

indir=$2

# SE ${indir} NÃO EXISTE, OU SEJA, SE NÃO FOI PASSADO ARGUMENTO 1 NA LINHA DE COMANDO
if [ ! ${indir} ]; then
	echo "[ERROR] Missing input directory." 1>&2
	exit
fi

# SE ${indir} NÃO É DIRETÓRIO
if [ ! -d ${indir} ]; then
	echo "[ERROR] Wrong input directory (${indir})." 1>&2
	exit
fi

outdir=$3

# SE ${outdir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 2 NA LINHA DE COMANDO
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
THREADS=$4

if [ ! ${THREADS} ]; then
	echo "[ERROR] Missing number of threads." 1>&2
	exit
fi

refgtf=$5
# SE ${refgtf} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 3 NA LINHA DE COMANDO
if [ ! ${refgtf} ]; then
	echo "[ERROR] Missing GTF file." 1>&2
	exit
fi

if [ ! -e "${refgtf}" ]; then
	echo "[ERROR] Not found GTF file (${refgtf})." 1>&2
	exit
fi

refseq=$6
# SE ${refseq} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 4 NA LINHA DE COMANDO
if [ ! ${refseq} ]; then
	echo "[ERROR] Missing GENOME fasta file." 1>&2
	exit
fi

if [ ! -e "${refseq}" ]; then
	echo "Not found GENOME fasta file (${refseq})." 1>&2
	exit
fi

# Opção cufflinks/stringtie
assembler=$7

if [ ! ${assembler} ]; then
	echo "[ERROR] Missing assembler (cufflinks or stringtie)." 1>&2
	exit
fi

if [ "${assembler}" != "cufflinks" ] &&
   [ "${assembler}" != "stringtie" ]; then

	echo "[ERROR] Assembler must be \"cufflinks\" or \"stringtie\" (${assembler})." 1>&2
  	exit
fi


# Contaminantes
contaminants=$8

./preprocess4.sh "${indir}" "${outdir}" "${THREADS}" ${contaminants}

echo -e "Starting Transcriptome Assembly ..."

# Criação de estrutura de diretórios
	
curdir=`pwd`
	
refseq_abs_path=$(readlink -f ${refseq})

if [ "${aligner}" == "tophat" ]; then

	mkdir -p ${outdir}/tophat_index
	mkdir -p ${outdir}/tophat_out_pe
	mkdir -p ${outdir}/tophat_out_se
	mkdir -p ${outdir}/tophat_out_final
	
	if [ ! -e "${outdir}/tophat_index/genome.fa" ]; then
		cd ${outdir}/tophat_index
		ln -s ${refseq_abs_path} genome.fa
		cd ${curdir}
	fi
		
	if [ ! -e "${outdir}/tophat_index/genome.1.bt2" ]; then
		echo -e "Indexing genome with TopHat2 ..."
		cd ${outdir}/tophat_index
		bowtie2-build 	--threads ${THREADS} \
				genome.fa genome > bowtie2.out.txt 2> bowtie2.err.txt
		cd ${curdir}
	fi
else
	# CASO CONTRÁRIO SERÁ star

	mkdir -p ${outdir}/star_index
	mkdir -p ${outdir}/star_out_pe
	mkdir -p ${outdir}/star_out_se
	mkdir -p ${outdir}/star_out_final
	
	if [ ! -e "${outdir}/star_index/genome.fa" ]; then
		cd ${outdir}/star_index
		ln -s ${refseq_abs_path} genome.fa
		cd ${curdir}
	fi
		
	if [ ! -e "${outdir}/star_index/SAindex" ]; then

		absrefgtf=`readlink -f ${refgtf}`

		cd ${outdir}/star_index
		
		echo -e "Indexing genome with STAR ..."
		
		# comando STAR para a indexação do genoma e transcriptoma
		STAR 	--runThreadN ${THREADS} \
	       		--runMode genomeGenerate \
			--genomeFastaFiles genome.fa \
			--genomeDir ./ \
			--sjdbGTFfile ${absrefgtf} \
			--genomeSAindexNbases 10 \
			--sjdbOverhang 149 \
			 > STAR.genomeGenerate.log.out.txt \
			2> STAR.genomeGenerate.log.err.txt

		cd ${curdir}
	fi

fi

if [ "${assembler}" == "cufflinks" ]; then
	mkdir -p ${outdir}/${aligner}_cufflinks
	mkdir -p ${outdir}/${aligner}_cuffmerge
elif [ "${assembler}" == "stringtie" ]; then
	mkdir -p ${outdir}/${aligner}_stringtie
	mkdir -p ${outdir}/${aligner}_stringmerge
fi

mkdir -p ${outdir}/${aligner}_${assembler}_cuffcompare
mkdir -p ${outdir}/${aligner}_${assembler}_cuffquant
mkdir -p ${outdir}/${aligner}_${assembler}_cuffnorm
mkdir -p ${outdir}/${aligner}_${assembler}_cuffdiff

for r1 in `find ${outdir}/ -name '*.prinseq.cleaned_1.fastq'`; do 

	r2=`echo ${r1} | sed 's/prinseq.cleaned_1.fastq/prinseq.cleaned_2.fastq/'`
	
	if [ ! -e ${r2} ]; then
		echo "[ERROR] Not found R2 (${r2})."
		exit
	fi
	
	echo -e "\tFound R1 ($(basename ${r1})) & R2 ($(basename ${r2})) ..."
	
	r1_singletons=`echo ${r1} | sed 's/prinseq.cleaned_1.fastq/prinseq.cleaned_1_singletons.fastq/'`
	r2_singletons=`echo ${r2} | sed 's/prinseq.cleaned_2.fastq/prinseq.cleaned_2_singletons.fastq/'`
	
	if [ ! -e ${r1_singletons} ]; then
		echo "[ERROR] Not found R1 singletons (${r1_singletons})."
		exit
	fi

	if   [ ! -e ${r2_singletons} ]; then
		echo "[ERROR] Not found R2 singletons (${r2_singletons})."
		exit
	fi
	
	name=`basename ${r1} .fastq | sed 's/.atropos_final.prinseq.cleaned_1//'`
	
	mkdir -p ${outdir}/align_out_final/${name}

	if [ "${aligner}" == "tophat" ]; then

		if [ ! -e "${outdir}/tophat_out_pe/${name}/accepted_hits.bam" ]; then

			echo -e "\tTopHat2 alignment (${name}) paired-end reads X genome ..." 

			tophat2 --min-anchor 75 \
				--min-intron-length 5000000 \
				--max-intron-length 5000000 \
				--no-novel-juncs \
				--max-multihits 20 \
				--transcriptome-max-hits 10 \
				--prefilter-multihits \
				--num-threads ${THREADS} \
				--GTF ${refgtf} \
				--transcriptome-index ${outdir}/tophat_index/transcriptome \
				--mate-inner-dist 0 \
				--mate-std-dev 100 \
				--coverage-search \
				--microexon-search \
				--b2-very-sensitive \
				--library-type fr-unstranded \
				--output-dir ${outdir}/tophat_out_pe/${name} \
				--no-sort-bam \
				${outdir}/tophat_index/genome \
				${r1} \
				${r2}	 > ${outdir}/tophat_out_pe/${name}.log.out.txt \
					2> ${outdir}/tophat_out_pe/${name}.log.err.txt
		else
			echo -e "\tFound Tophat2 output for PE (${name})..."
		fi		
		
		if [ ! -e "${outdir}/tophat_out_se/${name}/accepted_hits.bam" ]; then
			
			mkdir -p ${outdir}/tophat_out_se/${name}

			cat ${r1_singletons} ${r2_singletons} > ${outdir}/tophat_out_se/${name}/singletons.fastq
			
			if [ -s "${outdir}/tophat_out_se/${name}/singletons.fastq" ]; then

				echo -e "\tTopHat2 alignment (${name}) singleton reads X genome ..." 

				tophat2 --min-anchor 75 \
					--min-intron-length 5000000 \
					--max-intron-length 5000000 \
					--no-novel-juncs \
					--max-multihits 20 \
					--transcriptome-max-hits 10 \
					--prefilter-multihits \
					--num-threads ${THREADS} \
					--GTF ${refgtf} \
					--transcriptome-index ${outdir}/tophat_index/transcriptome \
					--coverage-search \
					--microexon-search \
					--b2-very-sensitive \
					--library-type fr-unstranded \
					--output-dir ${outdir}/tophat_out_se/${name} \
					--no-sort-bam \
					${outdir}/tophat_index/genome \
					${outdir}/tophat_out_se/${name}/singletons.fastq \
					  > ${outdir}/tophat_out_se/${name}.log.out.txt \
					 2> ${outdir}/tophat_out_se/${name}.log.err.txt

			fi
		else
			echo -e "\tFound Tophat2 output for SE (${name})..."
		fi
		
		if [ ! -e "${outdir}/tophat_out_final/${name}/accepted_hits.bam" ]; then
			
			mkdir -p ${outdir}/tophat_out_final/${name}
		
			if [ -s "${outdir}/tophat_out_pe/${name}/accepted_hits.bam" ]; then

				if [ -s "${outdir}/tophat_out_se/${name}/accepted_hits.bam" ]; then

					echo -e "\tMerging TopHat2 results ..."
		
					samtools view -H ${outdir}/tophat_out_pe/${name}/accepted_hits.bam > ${outdir}/tophat_out_final/${name}/Header.txt
					samtools merge 	-n --threads ${THREADS} \
							-h ${outdir}/tophat_out_final/${name}/Header.txt \
							${outdir}/tophat_out_final/${name}/accepted_hits.bam \
							${outdir}/tophat_out_pe/${name}/accepted_hits.bam \
							${outdir}/tophat_out_se/${name}/accepted_hits.bam \
						 > ${outdir}/tophat_out_final/${name}.log.out.txt \
						2> ${outdir}/tophat_out_final/${name}.log.err.txt

				else
					pe_result_abs_path=$(readlink -f ${outdir}/tophat_out_pe/${name}/accepted_hits.bam)
					cd ${outdir}/tophat_out_final/${name}/
					ln -s ${pe_result_abs_path} accepted_hits.bam
					cd ${curdir}
				fi				
			else
				if [ -s "${outdir}/tophat_out_se/${name}/accepted_hits.bam" ]; then
					se_result_abs_path=$(readlink -f ${outdir}/tophat_out_se/${name}/accepted_hits.bam)
					cd ${outdir}/tophat_out_final/${name}/
					ln -s ${se_result_abs_path} accepted_hits.bam
					cd ${curdir}
				else
					echo -e "[ERROR] Not found any alignment for PE or SE reads." 1>&2
				fi
			fi
		else
			echo -e "\tFound Tophat2 output final (${name})..."
		fi
					
		if [ ! -e "${outdir}/tophat_out_final/${name}/Aligned.sorted.bam" ]; then
			echo -e "\tSorting alignments (${name})..."
			samtools sort 	--threads ${THREADS} \
					-o ${outdir}/tophat_out_final/${name}/Aligned.sorted.bam \
			       		${outdir}/tophat_out_final/${name}/accepted_hits.bam \
					 > ${outdir}/tophat_out_final/${name}/Aligned.sorted.out.txt \
					2> ${outdir}/tophat_out_final/${name}/Aligned.sorted.err.txt
		fi
		
		# SEMPRE VAMOS REMOVER O LINK SIMBÓLICO PARA QUE AO ESCOLHER UM OUTRO 
		# ALINHADOR ELE SEJA SUBSTITUÍDO
		#if [ ! -e "${outdir}/align_out_final/${name}/Aligned.out.bam" ]; then
			rm -f ${outdir}/align_out_final/${name}/Aligned.out.bam
			rm -f ${outdir}/align_out_final/${name}/Aligned.sorted.bam
			if [ -e "${outdir}/tophat_out_final/${name}/accepted_hits.bam" ]; then
				align_final_out=`readlink -f ${outdir}/tophat_out_final/${name}/accepted_hits.bam`
				align_sorted_out=`readlink -f ${outdir}/tophat_out_final/${name}/Aligned.sorted.bam`
				cd ${outdir}/align_out_final/${name}
				ln -s ${align_final_out} Aligned.out.bam
				ln -s ${align_sorted_out} Aligned.sorted.bam
				cd ${curdir}
			else
				echo "[ERROR] Not found Tophat final output (${outdir}/tophat_out_final/${name}/accepted_hits.bam)" 2>&1
				exit
			fi				
		#fi
	else
		# SE NÃO FOR tophat ENTÃO star
		
		if [ ! -e "${outdir}/star_out_pe/${name}/Aligned.out.bam" ]; then

			echo -e "\tSTAR alignment (${name}) paired-end reads X genome ..." 

			mkdir -p ${outdir}/star_out_pe/${name}/

			STAR 	--runThreadN ${THREADS} \
				--genomeDir ${outdir}/star_index/ \
				--readFilesIn ${r1} ${r2} \
				--outSAMstrandField intronMotif \
				--outFilterIntronMotifs RemoveNoncanonical \
				--sjdbGTFfile ${refgtf} \
				--outFilterMultimapNmax 20 \
				--outFileNamePrefix ${outdir}/star_out_pe/${name}/ \
				--outSAMtype BAM Unsorted \
				--outFilterType BySJout \
				--outSJfilterReads Unique \
				--alignSJoverhangMin 75 \
				--alignSJDBoverhangMin 75 \
				--outFilterMismatchNmax 999 \
				--outFilterMismatchNoverReadLmax 0.04 \
				--alignIntronMin 5000000 \
				--alignIntronMax 1 \
				--alignMatesGapMax 1000000 \
			 > ${outdir}/star_out_pe/${name}.log.out.txt \
			2> ${outdir}/star_out_pe/${name}.log.err.txt
		else
			echo -e "\tFound STAR output for PE (${name})..."
		fi		

		if [ ! -e "${outdir}/star_out_se/${name}/Aligned.out.bam" ]; then
			
			mkdir -p ${outdir}/star_out_se/${name}

			cat ${r1_singletons} ${r2_singletons} > ${outdir}/star_out_se/${name}/singletons.fastq
			
			if [ -s "${outdir}/star_out_se/${name}/singletons.fastq" ]; then

				echo -e "\tSTAR alignment (${name}) singleton reads X genome ..." 
			
				STAR 	--runThreadN ${THREADS} \
					--genomeDir ${outdir}/star_index/ \
					--readFilesIn ${r1} \
					--outSAMstrandField intronMotif \
					--outFilterIntronMotifs RemoveNoncanonical \
					--sjdbGTFfile ${refgtf} \
					--outFilterMultimapNmax 20 \
					--outFileNamePrefix ${outdir}/star_out_se/${name}/ \
					--outSAMtype BAM Unsorted \
					--outFilterType BySJout \
					--outSJfilterReads Unique \
					--alignSJoverhangMin 75 
					--alignSJDBoverhangMin 75 \
					--outFilterMismatchNmax 999 \
					--outFilterMismatchNoverReadLmax 0.04 \
					--alignIntronMin 5000000 \
					--alignIntronMax 1 \
					--alignMatesGapMax 1000000 \
				 > ${outdir}/star_out_se/${name}.log.out.txt \
				2> ${outdir}/star_out_se/${name}.log.err.txt

			fi
		else
			echo -e "\tFound STAR output for SE (${name})..."
		fi
		
		if [ ! -e "${outdir}/star_out_final/${name}/Aligned.out.bam" ]; then
			
			mkdir -p ${outdir}/star_out_final/${name}
		
			if [ -s "${outdir}/star_out_pe/${name}/Aligned.out.bam" ]; then

				if [ -s "${outdir}/star_out_se/${name}/Aligned.out.bam" ]; then

					echo -e "\tMerging STAR results ..."
		
					samtools view -H ${outdir}/star_out_pe/${name}/Aligned.out.bam > ${outdir}/star_out_final/${name}/Header.txt
					samtools merge 	-n --threads ${THREADS} \
							-h ${outdir}/star_out_final/${name}/Header.txt \
							${outdir}/star_out_final/${name}/Aligned.out.bam \
							${outdir}/star_out_pe/${name}/Aligned.out.bam \
							${outdir}/star_out_se/${name}/Aligned.out.bam \
						 > ${outdir}/star_out_final/${name}.log.out.txt \
						2> ${outdir}/star_out_final/${name}.log.err.txt
		
		

				else
					pe_result_abs_path=$(readlink -f ${outdir}/star_out_pe/${name}/Aligned.out.bam)
					cd ${outdir}/star_out_final/${name}/
					ln -s ${pe_result_abs_path} Aligned.out.bam
					cd ${curdir}
				fi
			else
				if [ -s "${outdir}/star_out_se/${name}/Aligned.out.bam" ]; then
					se_result_abs_path=$(readlink -f ${outdir}/star_out_se/${name}/Aligned.out.bam)
					cd ${outdir}/star_out_final/${name}/
					ln -s ${se_result_abs_path} Aligned.out.bam
					cd ${curdir}
				else
					echo -e "[ERROR] Not found any alignment for PE or SE reads." 1>&2
				fi
			fi
		else
			echo -e "\tFound STAR output final (${name})..."
		fi
					
		if [ ! -e "${outdir}/star_out_final/${name}/Aligned.sorted.bam" ]; then

			echo -e "\tSorting alignments (${name})..."

			samtools sort 	--threads ${THREADS} \
					-o ${outdir}/star_out_final/${name}/Aligned.sorted.bam \
			       		${outdir}/star_out_final/${name}/Aligned.out.bam \
					 > ${outdir}/star_out_final/${name}/Aligned.sorted.out.txt \
					2> ${outdir}/star_out_final/${name}/Aligned.sorted.err.txt
		fi
		
		# SEMPRE VAMOS REMOVER O LINK SIMBÓLICO PARA QUE AO ESCOLHER UM OUTRO 
		# ALINHADOR ELE SEJA SUBSTITUÍDO
		#if [ ! -e "${outdir}/align_out_final/${name}/Aligned.out.bam" ]; then
			rm -f ${outdir}/align_out_final/${name}/Aligned.out.bam
			rm -f ${outdir}/align_out_final/${name}/Aligned.sorted.bam
			if [ -e "${outdir}/star_out_final/${name}/Aligned.out.bam" ]; then
				align_final_out=`readlink -f ${outdir}/star_out_final/${name}/Aligned.out.bam`
				align_sorted_out=`readlink -f ${outdir}/star_out_final/${name}/Aligned.sorted.bam`
				cd ${outdir}/align_out_final/${name}
				ln -s ${align_final_out} Aligned.out.bam
				ln -s ${align_sorted_out} Aligned.sorted.bam
				cd ${curdir}
			else
				echo "[ERROR] Not found STAR final output (${outdir}/star_out_final/${name}/Aligned.out.bam)" 2>&1
				exit
			fi
		#fi

	fi
	
	mkdir -p ${outdir}/align_out_info/
	
	if [ -e "${outdir}/align_out_final/${name}/Aligned.out.bam" ]; then
	
		if [ ! -e "${outdir}/align_out_info/${name}.${aligner}.out.txt" ]; then

			echo -e "\tGet alignment information (${name}) [${aligner}] ..."
			
			SAM_nameSorted_to_uniq_count_stats.pl ${outdir}/align_out_final/${name}/Aligned.out.bam > ${outdir}/align_out_info/${name}.${aligner}.out.txt 2> ${outdir}/align_out_info/${name}.${aligner}.err.txt
			
		fi
			
		mkdir -p ${outdir}/${aligner}_cufflinks/${name}

		if [ "${assembler}" == "cufflinks" ]; then
			if [ ! -e "${outdir}/${aligner}_cufflinks/${name}/transcripts.gtf" ]; then		
				echo -e "\tAssembly transcriptome (${name}) [cufflinks]\n"

				cufflinks	--output-dir ${outdir}/${aligner}_cufflinks/${name} \
						--num-threads ${THREADS} \
						--GTF ${refgtf} \
						--frag-bias-correct ${refseq} \
						--multi-read-correct \
						--library-type fr-unstranded \
						--frag-len-mean 300 \
						--frag-len-std-dev 30 \
						--total-hits-norm \
						--min-isoform-fraction 0.5 \
						--pre-mrna-fraction 0.5 \
						--min-frags-per-transfrag 10 \
						--junc-alpha 0.0000009 \
						--small-anchor-fraction 0.5 \
						--overhang-tolerance 0 \
						--min-intron-length 5000000 \
						--max-intron-length 1 \
						--trim-3-avgcov-thresh 0.05 \
						--trim-3-dropoff-frac 0.01 \
						--max-multiread-fraction 0.75 \
						--overlap-radius 1 \
						--3-overhang-tolerance 0 \
						--intron-overhang-tolerance 0 \
						${outdir}/align_out_final/${name}/Aligned.sorted.bam \
					 > ${outdir}/${aligner}_cufflinks/${name}/cufflinks.out.txt \
					2> ${outdir}/${aligner}_cufflinks/${name}/cufflinks.err.txt
			fi

		else
			if [ ! -e "${outdir}/${aligner}_stringtie/${name}/transcripts.gtf" ]; then		
				mkdir -p ${outdir}/${aligner}_stringtie/${name}
					
				echo -e "\tAssembly transcriptome (${name}) [stringtie]\n"

				stringtie	${outdir}/align_out_final/${name}/Aligned.sorted.bam \
						-G ${refgtf} \
						-f 0.5 \
						-m 200 \
						-o ${outdir}/${aligner}_stringtie/${name}/transcripts.gtf \
						-a 75 \
						-j 10000 \
						-c 4 \
						-v \
						-g 0 \
						-C ${outdir}/${aligner}_stringtie/${name}/coverages.txt \
						-M 1 \
						-p ${THREADS} \
						-A ${outdir}/${aligner}_stringtie/${name}/abundances.txt \
						-B \
					 > ${outdir}/${aligner}_stringtie/${name}/stringtie.out.txt \
					2> ${outdir}/${aligner}_stringtie/${name}/stringtie.err.txt				
			fi

		fi


	else
		echo -e "[ERROR] Not found alignment data (${outdir}/align_out_final/${name}/Aligned.out.bam)"	1>&2
	fi
done

rm -f ${outdir}/assembly_GTF_list.txt
for transc in `find ${outdir}/${aligner}_${assembler} -name transcripts.gtf`; do
	#echo -e "\tProcessing transcriptome ${transc} ..."
	echo ${transc} >> ${outdir}/assembly_GTF_list.txt
done

transcriptomeref=""

if [ "${assembler}" == "cufflinks" ]; then

	if [ ! -e "${outdir}/${aligner}_cuffmerge/merged.gtf" ]; then

		echo -e "\tMerging transcriptomes (${outdir}/assembly_GTF_list.txt) in a transcriptome reference [cuffmerge]\n"
	
		cuffmerge 	-o ${outdir}/${aligner}_cuffmerge \
				--ref-gtf ${refgtf} \
				--ref-sequence ${refseq} \
				--min-isoform-fraction 0.5 \
				--num-threads ${THREADS} \
				${outdir}/assembly_GTF_list.txt \
			 > ${outdir}/${aligner}_cuffmerge/cuffmerge.out.txt \
			2> ${outdir}/${aligner}_cuffmerge/cuffmerge.err.txt
	fi
	
	transcriptomeref="${outdir}/${aligner}_cuffmerge/merged.gtf"
else

	if [ ! -e "${outdir}/${aligner}_stringmerge/merged.gtf" ]; then

		echo -e "\tMerging transcriptomes (${outdir}/assembly_GTF_list.txt) in a transcriptome reference [stringtie]\n"
		
		stringtie 	--merge	\
				-G ${refgtf} \
				-o ${outdir}/${aligner}_stringmerge/merged.gtf \
				-m 200 \
				-c 1 \
				-F 4 \
				-T 4 \
				-f 0.5 \
				-g 0 \
				${outdir}/assembly_GTF_list.txt \
			> ${outdir}/${aligner}_stringmerge/stringmerge.out.txt \
			2> ${outdir}/${aligner}_stringmerge/stringmerge.err.txt

	fi
		
	transcriptomeref="${outdir}/${aligner}_stringmerge/merged.gtf"
	
	
fi
	
if [ ! -e "${outdir}/${aligner}_${assembler}_cuffcompare/cuffcmp.combined.gtf" ]; then
	
	echo -e "\tRunning cuffcompare with ${aligner} & ${assembler} transcriptome reference (${transcriptomeref})..."

	cuffcompare 	-r ${refgtf} \
			-s ${refseq} \
			-o ${outdir}/${aligner}_${assembler}_cuffcompare/cuffcmp \
			${transcriptomeref} \
		 > ${outdir}/${aligner}_${assembler}_cuffcompare/cuffcmp.out.txt \
		2> ${outdir}/${aligner}_${assembler}_cuffcompare/cuffcmp.err.txt

fi

# LISTA DE VALORES NÃO REDUNDANTES (NOME DO GRUPO BIOLÓGICO)
# Ex.: (CONTROL TEST)
biogroup_label=()

for bamfile in `find ${outdir}/align_out_final -name Aligned.sorted.bam`; do

	name=`basename $(dirname ${bamfile})`

	if [ ! -e "${outdir}/${aligner}_${assembler}_cuffquant/${name}/abundances.cxb" ]; then

		echo -e "\tRunning cuffquant using sample ${name} as using ${aligner} & ${assembler} (${transcriptomeref}) ..."
		mkdir -p ${outdir}/${aligner}_${assembler}_cuffquant/${name}
			
		cuffquant 	--output-dir ${outdir}/${aligner}_${assembler}_cuffquant/${name} \
				--frag-bias-correct ${refseq} \
				--multi-read-correct \
				--num-threads ${THREADS} \
				--library-type fr-unstranded \
				--frag-len-mean 300 \
				--frag-len-std-dev 30 \
				--max-bundle-frags 9999999 \
				--max-frag-multihits 20 \
				${transcriptomeref} \
				${bamfile} \
			 > ${outdir}/${aligner}_${assembler}_cuffquant/${name}/cuffquant.log.out.txt \
			2> ${outdir}/${aligner}_${assembler}_cuffquant/${name}/cuffquant.log.err.txt

	fi
	
	groupname=`echo ${name} | sed 's/[0-9]\+$//'`
	biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))

done

biogroup_files=()

echo -e "\tCollecting Expression Data from cuffquant output (*.cxb) ..."

for label in ${biogroup_label[@]}; do
	echo -e "\t\tCollecting .cxb files for ${label} ..."
	group=()
	for cxbfile in `ls ${outdir}/${aligner}_${assembler}_cuffquant/${label}*/abundances.cxb`; do
		echo -e "\t\t\tFound ${cxbfile}"
		group=(${group[@]} "${cxbfile}")
	done
	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "Starting Gene Expression Analysis ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}

if [ ! -e "${outdir}/${aligner}_${assembler}_cuffnorm/isoforms.count_table" ]; then

	echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

	cuffnorm 	--output-dir ${outdir}/${aligner}_${assembler}_cuffnorm \
			--labels $(IFS=, ; echo "${biogroup_label[*]}") \
			--num-threads ${THREADS} \
			--library-type fr-unstranded \
			--library-norm-method geometric \
			--output-format simple-table \
			${transcriptomeref} \
			${biogroup_files[*]} \
			 > ${outdir}/${aligner}_${assembler}_cuffnorm/cuffnorm.log.out.txt \
			2> ${outdir}/${aligner}_${assembler}_cuffnorm/cuffnorm.log.err.txt

fi

if [ ! -e "${outdir}/${aligner}_${assembler}_cuffnorm/isoforms.raw_count_table.txt" ]; then
	de-normalize-cuffnorm.R --in=${outdir}/${aligner}_${assembler}_cuffnorm/isoforms.count_table \
				--st=${outdir}/${aligner}_${assembler}_cuffnorm/samples.table \
				--out=${outdir}/${aligner}_${assembler}_cuffnorm/isoforms.raw_count_table.txt
	 > ${outdir}/${aligner}_${assembler}_cuffnorm/de-normalize-cuffnorm.out.txt \
	2> ${outdir}/${aligner}_${assembler}_cuffnorm/de-normalize-cuffnorm.err.txt

fi

if [ ! -e "${outdir}/${aligner}_${assembler}_cuffdiff/isoform_exp.diff" ]; then
	echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

	cuffdiff 	--output-dir ${outdir}/${aligner}_${assembler}_cuffdiff \
			--labels $(IFS=, ; echo "${biogroup_label[*]}") \
			--frag-bias-correct ${refseq} \
			--multi-read-correct \
			--num-threads ${THREADS} \
			--library-type fr-unstranded \
			--frag-len-mean 300 \
			--frag-len-std-dev 30 \
			--max-bundle-frags 9999999 \
			--max-frag-multihits 20 \
			--total-hits-norm \
			--min-reps-for-js-test 2 \
			--library-norm-method geometric \
			--dispersion-method per-condition \
			--min-alignment-count 10 \
			${transcriptomeref} \
			${biogroup_files[*]} \
			 > ${outdir}/${aligner}_${assembler}_cuffdiff/cuffdiff.log.out.txt \
			2> ${outdir}/${aligner}_${assembler}_cuffdiff/cuffdiff.log.err.txt
fi

